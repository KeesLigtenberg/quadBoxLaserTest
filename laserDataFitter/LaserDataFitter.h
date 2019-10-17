/*
 * LaserDataFitter.h
 *
 *  Created on: Jul 16, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_LASERDATAFITTER_H_
#define LASERDATAFITTER_LASERDATAFITTER_H_

#include <string>
#include <vector>
#include <memory>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

#include "../eventBuilder/Hit.h"
#include "PositionHit.h"
#include "Alignment.h"

struct QuadTreeReader {
	QuadTreeReader(std::string treeName, std::string fileName) :
				file(new TFile(fileName.c_str(), "READ")),
				tree(nullptr) {
		tree= dynamic_cast<TTree*>(file->Get(treeName.c_str()));
		if(!tree) {
			std::cerr<<"could not get tree "<<treeName<<" from file "<<file->GetName()<<std::endl;
			throw 1;
		}

		tree->SetBranchAddress("triggerToA", &triggerToA);
		tree->SetBranchAddress("triggerNumber", &triggerNumber);
		tree->SetBranchAddress( "laser", &laser );
		for(int i=0; i<4; i++) {
			tree->SetBranchAddress( ("chip"+std::to_string(i)).c_str(), chip+i );
		}
	}

	Long64_t entryNumber{0};
	bool getNextEntry(int increaseEntryBy=1) {
		entryNumber+=increaseEntryBy; //change this for quick alignment
		static Long64_t nEntries=tree->GetEntries();
		if(!(entryNumber%1000)) std::cout<<entryNumber<<"/"<<nEntries<<"\n";
		if(entryNumber < nEntries) {
			return tree->GetEntry(entryNumber);
		} else return false;
	}

	std::unique_ptr<TFile> file;
	TTree* tree;

	Long64_t triggerToA{};
  unsigned triggerNumber{};
  Vec3* laser{};
	std::vector<Hit>* chip[4]{};
};



class QuadTrackFitter {
public:
	QuadTrackFitter(std::string fileName);
	virtual ~QuadTrackFitter();
	void Loop(std::string outputFile, const Alignment& alignment);

	QuadTreeReader reader;
	std::vector<PositionHit> posHits;
	Vec3 laser;

	std::function<bool(const Vec3&)> selectLaserPoint=[](const Vec3&){return true;};

};

#endif /* LASERDATAFITTER_LASERDATAFITTER_H_ */
