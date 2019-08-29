/*
 * detectorConfiguration.h
 *
 *  Created on: Jun 23, 2017
 *      Author: cligtenb
 */

#ifndef DETECTORCONFIGURATION_H_
#define DETECTORCONFIGURATION_H_

#include <utility>
#include "TVector3.h"

struct DetectorConfiguration {
	virtual ~DetectorConfiguration() {};

	virtual double xmax() const =0;
	virtual double xmin() const { return 0.; }
	virtual double ymax() const  =0;
	virtual double ymin() const { return 0.; }
	virtual double zmax() const  =0;
	virtual double zmin() const { return 0; }
	virtual double zmean() const { return (zmin()+zmax())/2; }

};

struct SimpleDetectorConfiguration : DetectorConfiguration {
	SimpleDetectorConfiguration(double minx, double maxx, double miny, double maxy, double minz, double maxz) : minx(minx), maxx(maxx), miny(miny), maxy(maxy), minz(minz), maxz(maxz) {};
	double minx, maxx, miny, maxy, minz, maxz;
	virtual ~SimpleDetectorConfiguration() {};
	virtual double xmax() const { return maxx; }
	virtual double xmin() const { return minx; }
	virtual double ymax() const { return maxy; }
	virtual double ymin() const { return miny; }
	virtual double zmax() const { return maxz; }
	virtual double zmin() const { return minz; }
};

SimpleDetectorConfiguration simpleDetectorFromChipCorners( std::vector<TVector3> corners) {
	auto xmm = std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w) {return v.x()<w.x();});
	auto ymm = std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w) {return v.y()<w.y();});
	auto zmm = std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w) {return v.z()<w.z();});
	return {xmm.first->x(), xmm.second->x(), ymm.first->y(), ymm.second->y(),zmm.first->z(), zmm.second->z()};
}



#endif /* DETECTORCONFIGURATION_H_ */
