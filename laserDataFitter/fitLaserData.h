/*
 * fitLaserData.h
 *
 *  Created on: Jul 26, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_FITLASERDATA_H_
#define LASERDATAFITTER_FITLASERDATA_H_

#include "LaserDataFitter.cpp"
#include "/user/cligtenb/rootmacros/AllCombiner.h"


struct AlignmentGrapher {

	TGraph positionx, positiony, positionz;
	TGraph rotationx, rotationy, rotationz;
	int iteration=0;
	TVector3 firstShift{}, firstRotation{};

	void addPoint(const Alignment& a) {
		if(not iteration) {
				firstShift=a.quad.shift;
				firstRotation=a.quad.rotation;
		}
		positionx.SetPoint(iteration, iteration, a.quad.shift.x()-firstShift.x());
		positiony.SetPoint(iteration, iteration, a.quad.shift.y()-firstShift.y());
		positionz.SetPoint(iteration, iteration, a.quad.shift.z()-firstShift.z());
		rotationx.SetPoint(iteration, iteration, a.quad.rotation.x()-firstRotation.x());
		rotationy.SetPoint(iteration, iteration, a.quad.rotation.y()-firstRotation.y());
		rotationz.SetPoint(iteration, iteration, a.quad.rotation.z()-firstRotation.z());
		iteration++;
	}

	TCanvas* createCanvas() {
		AllCombiner<TGraph> position{"combinedPosition"};
		position.add((TGraph*)positionx.Clone(), "Quad x-position");
		position.add((TGraph*)positiony.Clone(), "Quad y-position");
		position.add((TGraph*)positionz.Clone(), "Quad z-position");
		position.setStyle(7);

		AllCombiner<TGraph> rotation{"combinedRotation"};
		rotation.add((TGraph*)rotationx.Clone(), "Quad x-rotation");
		rotation.add((TGraph*)rotationy.Clone(), "Quad y-rotation");
		rotation.add((TGraph*)rotationz.Clone(), "Quad z-rotation");
		rotation.setStyle(7);

		rotation.createCombined();
		return position.createCombined();
	}

};


void fitLaserData(std::string inputFileName, std::string outputFileName, std::string alignFile="align.dat") {
	QuadTrackFitter ldf{inputFileName};
	Alignment alignment{alignFile};
	//do fit
	ldf.Loop(outputFileName, alignment);
}

void fitAndUpdateAlignment(std::string inputFileName, std::string outputFileName, std::string alignFile="align.dat", int nRepeat=1) {

	AlignmentGrapher alignmentGraphs;

	for(int i=0; i<nRepeat; i++)  {
		std::cout<<"iteration "<<i<<"\n\n";

		QuadTrackFitter ldf{inputFileName};
//		ldf.minDistanceFromEdge=1.3;//mm
//		ldf.selectLaserPoint=minZLaserCut{9}; //9 mm
//		ldf.selectLaserPoint=cutByMinimumDistanceFromEdge{1.3};
		Alignment alignment{alignFile};

		alignmentGraphs.addPoint(alignment);

		//do fit
		ldf.Loop(outputFileName, alignment);

		TFile file(outputFileName.c_str(), "READ");
//		alignment.updateAll(file);
		alignment.updateShifts(file);
//		alignment.quad.updateShift(file, "quad/global");
		alignment.updateRotations(file);
//		alignment.quad.updateRotation(file, "quad");
//		alignment.timeWalk.update(file);
		alignment.updateDriftSpeed(file);
		alignment.write(alignFile);

	}

	alignmentGraphs.createCanvas();

}


#endif /* LASERDATAFITTER_FITLASERDATA_H_ */
