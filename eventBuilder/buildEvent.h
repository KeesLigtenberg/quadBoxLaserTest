/*
d * buildEvent.h
 *
 *  Created on: Oct 8, 2018
 *      Author: cligtenb
 */

#ifndef EVENTBUILDER_BUILDEVENT_H_
#define EVENTBUILDER_BUILDEVENT_H_


#include <string>
#include <array>
#include <vector>
#include <memory>
#include <algorithm>
#include <deque>
#include <list>

#include "TTree.h"
#include "TH1.h"
#include "TProfile2D.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"

//#include "TriggerDecoder.h" for testbeam
#include "Hit.h"


//run name from folder name
std::string getRunFromFolder(std::string runDir) {
	auto pos=runDir.rfind('/');
	if(pos!=std::string::npos) {//found
		return runDir.substr(pos);
	} else {
		return runDir;
	}
}

struct TreeReader {
	TTree* tree{};
	unsigned currentEntry{0};
	const unsigned nEntries{0};

	TreeReader() {};
	TreeReader(TTree* tree) : tree(tree), nEntries(tree->GetEntries()) {};

	virtual ~TreeReader() {};

	virtual void getFirstEntry() {
		currentEntry=0;
		tree->GetEntry(currentEntry);
	}

	virtual void getNextEntry() {
		tree->GetEntry(++currentEntry);
	}
	bool reachedEnd() {
		return currentEntry+1>=nEntries;
	}
};


struct ChipTreeReader : TreeReader {
	unsigned long long toa{};
	unsigned char col{},row{};
	unsigned short tot{};

	ChipTreeReader( TTree* tree ) : TreeReader(tree) {
		tree->SetBranchAddress("toa", &toa);
		tree->SetBranchAddress("col", &col);
		tree->SetBranchAddress("row", &row);
		tree->SetBranchAddress("tot", &tot);
		tree->GetEntry(currentEntry);
	}
	virtual void getNextEntry() {
		auto oldToa=toa;
		TreeReader::getNextEntry();
		if(toa<oldToa) std::cout<<"Warning: toa decreased by "<<oldToa-toa<<" to "<<toa<<"\n";
	}
	void discardEntriesBeforeT(unsigned long long t) {
		while(not reachedEnd() and toa < t) getNextEntry();
	}
};

struct TriggerTreeReader : TreeReader {
	unsigned long long toa{};
	TriggerTreeReader( TTree* tree ) : TreeReader(tree) {
		tree->SetBranchAddress("timestamp", &toa);
		tree->GetEntry(currentEntry);
	}
};

struct StageTreeReader : TreeReader {
	long long toa{};
	double pos[3]{};

	StageTreeReader( TTree* tree ) : TreeReader(tree) {
		tree->SetBranchAddress("toa", &toa);
		tree->SetBranchAddress("x", pos);
		tree->SetBranchAddress("y", pos+1);
		tree->SetBranchAddress("z", pos+2);
	}

	~StageTreeReader() {
		delete tree;
	}

	virtual void getNextEntry() {
		auto previousToa=toa;
		tree->GetEntry(++currentEntry);

		if(not toa) { //toa was not read correctly
			tree->GetEntry(currentEntry+1);
			auto nextToa=toa;
			if(not nextToa) {
				std::cerr<<"two stage positions were missing in sequence!";
				throw;
			}
			tree->GetEntry(currentEntry);
			toa=(nextToa+previousToa)/2; //place toa between next and previous
			std::cout<<"missing toa, put it at "<<toa<<" between "<<nextToa<<" and "<<previousToa<<"\n";
		}
	}

};


struct CombinedTreeWriter {
	TTree tree;
	std::vector<Hit> chips[4]{};
	Vec3	laser; //laserposition
	long long triggerToA=0;
	unsigned triggerNumber=0;

	CombinedTreeWriter() : tree{"data", "tree with run data"} {
		tree.Branch( "triggerToA", &triggerToA );
		tree.Branch( "laser", &laser );
		tree.Branch( "triggerNumber", &triggerNumber );
		for(int i=0; i<4; i++)
			tree.Branch( ("chip"+std::to_string(i)).c_str() , "std::vector<Hit>", chips+i);
	}
	void fill() {
			tree.Fill();
	}
	void clear() {
		for(auto& c : chips) c.clear();
	}
};

void convertToTree(std::string inputFileName, std::string stagePositionFileName, std::string outputFileName) {

	//Read trees
	std::cout<<"reading trees\n";
	std::vector<std::unique_ptr<ChipTreeReader> > chips;
	TFile inputFile(inputFileName.c_str(), "READ");
	for(int i=0; i<4; i++) {
		std::cout<<"  chip "<<i<<"\n";
		auto chipTree=getObjectFromFile<TTree>("hits_chip"+std::to_string(i), &inputFile);
		chips.emplace_back( new ChipTreeReader(chipTree) );
	}
	std::cout<<"  trigger\n";
	auto triggerTree=getObjectFromFile<TTree>("triggers", &inputFile);
	TriggerTreeReader trigger(triggerTree);
	std::cout<<"  stage\n";
	auto stageTree=new TTree();
	stageTree->ReadFile(stagePositionFileName.c_str(), "toa/L:x/D:y:z");
	StageTreeReader stage(stageTree);
	stage.getFirstEntry();

	//Output tree
	TFile file((outputFileName).c_str(), "RECREATE"); //use absolute path
	CombinedTreeWriter output;

	bool stagePositionSet=false;
	while( not trigger.reachedEnd() ) {
		//find next signal
		if(!(trigger.currentEntry%1000) ) {
			std::cout<<"trigger "<<trigger.currentEntry<<"/ "<<trigger.tree->GetEntries()<<": "<<trigger.toa<<", stage "<<stage.currentEntry<<": "<<stage.toa*4096<<"\n";
		}

		if(Long64_t(trigger.toa) > stage.toa*4096 and not stage.reachedEnd() ) { //stage first
			//move to next stage position
			output.laser=Vec3(stage.pos[0], stage.pos[1], stage.pos[2]);
			stagePositionSet=true;
			stage.getNextEntry();
		} else { //trigger first
			//if no stage position yet, continue
			if(not stagePositionSet) {
				continue;
			}

			//set trigger
			output.triggerToA=trigger.toa;

			//find all hits within a range of the trigger
			const double maxTimeBeforeTrigger=500*4096/25; //500 ns
			const double maxTimeAfterTrigger=500*4096/25; //500 ns
			for(unsigned i=0; i<chips.size(); i++) {
				auto& c = chips[i];
				c->discardEntriesBeforeT(trigger.toa-maxTimeBeforeTrigger);
				while(c->toa < trigger.toa+maxTimeAfterTrigger and not c->reachedEnd()) {
					output.chips[i].emplace_back(c->row, c->col, c->tot, c->toa-trigger.toa );
					c->getNextEntry();
				}
			}//for chips

			//fill tree after each trigger
			output.fill();
			output.clear();
			trigger.getNextEntry();
		}
	} //while triggers and stage

	output.tree.Write("data");
	file.Write();
}


void buildEvent() {

}




#endif /* EVENTBUILDER_BUILDEVENT_H_ */
