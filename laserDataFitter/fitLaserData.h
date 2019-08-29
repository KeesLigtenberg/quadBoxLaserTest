/*
 * fitLaserData.h
 *
 *  Created on: Jul 26, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_FITLASERDATA_H_
#define LASERDATAFITTER_FITLASERDATA_H_

#include "LaserDataFitter.cpp"

void fitLaserData(std::string inputFileName, std::string outputFileName, std::string alignFile="align.dat") {
	QuadTrackFitter ldf{inputFileName};
	Alignment alignment{alignFile};
	//do fit
	ldf.Loop(outputFileName, alignment);
}

void fitAndUpdateAlignment(std::string inputFileName, std::string outputFileName, std::string alignFile="align.dat", int nRepeat=1) {

	for(int i=0; i<nRepeat; i++)  {
		std::cout<<"iteration "<<i<<"\n\n";

		QuadTrackFitter ldf{inputFileName};
//		ldf.minDistanceFromEdge=1.3;//mm
		ldf.selectLaserPoint=minZLaserCut{9}; //9 mm
		Alignment alignment{alignFile};

		//do fit
		ldf.Loop(outputFileName, alignment);

		TFile file(outputFileName.c_str(), "READ");
//		alignment.updateAll(file);
//		alignment.updateShifts(file);
		alignment.quad.updateShift(file, "quad/global");
//		alignment.updateRotations(file);
		alignment.quad.updateRotation(file, "quad");
//		alignment.timeWalk.update(file);
		alignment.updateDriftSpeed(file);
		alignment.write(alignFile);

	}

}


#endif /* LASERDATAFITTER_FITLASERDATA_H_ */
