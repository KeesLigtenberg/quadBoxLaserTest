/*
 * LaserDataFitter.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: cligtenb
 */

#include "LaserDataFitter.h"

#include <algorithm>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

#include "PositionHit.h"
#include "DetectorConfiguration.h"
#include "drawFunctions.h"

#pragma link C++ class std::vector<int>+;

namespace {

	const int nChips=4;
//	struct shift {double x,y;};
//	std::array<shift,nChips> chipShifts={{ {24.1,15.7}, {38.0, 26.3}, {23.7,26.3}, {9.9,15.7} }};


	//selection functions for laser points
	bool selectByRectangularCut(const Vec3& v) {
		struct {	double x,y;	} chipCenters[4] = { {17,8.5}, {31,8.5}, {17,33.5}, {31,33.5} } ;
		double maxDistance=5.1;
		for(auto& center : chipCenters) {
			if(std::fabs(center.x-v.x)<maxDistance and std::fabs(center.y-v.y)<maxDistance) {
				return true;
			}
		}
		return false;
	}
	struct cutByMinimumDistanceFromEdge{
		const Alignment& alignment;
		double minDistance;//mm
		bool operator() (const Vec3& laser) {
			if(std::any_of(alignment.chips.begin(), alignment.chips.end(), [&laser, this](const ChipAlignment& ca){
					return ca.getDistanceFromEdge(laser)<minDistance;
				})) {
				return false;
			}
			return true;
		}
	};
	struct minZLaserCut {
		double minZ;
		bool operator() (const Vec3& laser) {
			return laser.z>minZ;
		}
	};

	std::vector<PositionHit> convertHitsQuad( std::vector<Hit>* chips [], double driftSpeed, double t0Offset=0 /*value to subtract from time*/) {
		std::vector<PositionHit> hits;
		for(int i=0; i<nChips; i++) {
			auto chipHits=convertHitsTPC(*chips[i], i, driftSpeed, t0Offset); //t0 is subtracted from hit time
			hits.insert(hits.end(), chipHits.begin(), chipHits.end() );
		}
		return hits;
	}

	struct PositionRange { double xmin, xmax, ymin, ymax; };
	inline PositionRange rangeFromCorners(const std::vector<TVector3>& corners) {
		auto xrange=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w){return v.x()<w.x();});
		auto yrange=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w){return v.y()<w.y();});
		auto xmin=xrange.first->x(),xmax=xrange.second->x(),ymin=yrange.first->y(),ymax=yrange.second->y();
		return {xmin,xmax,ymin,ymax};
	}

}//namespace




struct ChipHistogrammer {
	ChipHistogrammer(const char* name, const Alignment& align) : ChipHistogrammer(std::string(name), align) {}
	ChipHistogrammer(std::string name, const Alignment& align) ;

	std::unique_ptr<TH1I> nHits;
	std::unique_ptr<TH1D> driftTime, ToT;
	std::unique_ptr<TH1D> xResidualPull{}, yResidualPull{}, zResidualPull{};
	std::unique_ptr<TH2D> pixelHitMap;
	std::unique_ptr<TH1D> xRotation, yRotation, zRotation;
	std::unique_ptr<TH1D> xzShear; //angle of E field to plane

	std::unique_ptr<TH2D> zResidualByToT, xResidualByToT, zResidualByToTCorrected,zResidualByDriftTime;
	std::unique_ptr<TH2D> xByzForShear;

	struct residualHistograms {
		residualHistograms(std::string name) : name{name} {};
		std::string name;
		std::unique_ptr<TH1D> xResidual{}, yResidual{}, zResidual{};
		std::unique_ptr<TH1D> zHit{};
		std::unique_ptr<TH2D> xResidualByz{}, yResidualByz{},zResidualByz{};
		std::unique_ptr<TH3D> zResidualByzByToT{};
		std::unique_ptr<TH2D>	positionHitMap{}, XZHitMap{},YZHitMap{};
		std::unique_ptr<TProfile2D> xResidualByPosition{}, zResidualByPosition{};
		std::unique_ptr<TProfile2D> xResidualByXZ{}, zResidualByXZ{};
		void fill(const Vec3& position, const Vec3& residual, double ToT) {
			xResidual->Fill(residual.x);
			yResidual->Fill(residual.y);
			zResidual->Fill(residual.z);
			xResidualByz->Fill(position.z,residual.x);
			yResidualByz->Fill(position.z,residual.y);
			zResidualByz->Fill(position.z,residual.z);
			zResidualByzByToT->Fill(position.z,residual.z, ToT*.025);
			positionHitMap->Fill(position.x, position.y);
			XZHitMap->Fill(position.x,position.z);
			YZHitMap->Fill(position.y,position.z);
			zHit->Fill(position.z);
			xResidualByPosition->Fill(position.x,position.y,residual.x);
			zResidualByPosition->Fill(position.x,position.y,residual.z);

			xResidualByXZ->Fill(position.x,position.z,residual.x);
			zResidualByXZ->Fill(position.x,position.z,residual.z);
		}
		void createHistograms(const PositionRange& range);
	} global{"_global"}, local{"_local"}, locExp{"_locExp"};
	std::unique_ptr<TProfile2D> xResidualByPixel, zResidualByPixel;

	void fillHit(const PositionHit& h);
	void fillRotation(const TVector3& pos, const TVector3& residual, const TVector3& COM);
	void fillRotation(const PositionHit& h, const TVector3& COM) { fillRotation(h.position, h.residual, COM); }
	void fillEvent();
	void fillTimewalk(double localZResidual, double ToT, const TimeWalkCorrector& twc) {
		zResidualByToT->Fill(ToT*25E-3, localZResidual+twc.getCorrection(ToT));
	}
private:
	std::string dirName;
	int hitsCounter=0;
};

ChipHistogrammer::ChipHistogrammer(std::string name, const Alignment& align) : dirName(name) {
	using namespace std;

	//change directory
	auto startDir=gDirectory;
	gDirectory->mkdir( dirName.c_str() )->cd();

	nHits=		unique_ptr<TH1I>(new TH1I{"nHits", "Number of hits per event;Number of hits;Entries", 200,0,400});
	constexpr double ToABinWidth=1.5625E-3;
	driftTime=unique_ptr<TH1D>(new TH1D{"driftTime","Drift time;Drift time [#mus];Hits", int(0.8/ToABinWidth),-0.39999,0.4}); //zero would be rounded up (wrong), so move bin slightly to 0.099999
	ToT=			unique_ptr<TH1D>(new TH1D{"ToT", "Time over threshold;ToT [#mus];Hits",80,0,2});

	xResidualPull=unique_ptr<TH1D>(new TH1D{("xResidualPull"), "Pull of x-residual;x-residual pull; Hits", 100,-4,4});
	yResidualPull=unique_ptr<TH1D>(new TH1D{("yResidualPull"), "Pull of y-residual;y-residual pull; Hits", 100,-4,4});
	zResidualPull=unique_ptr<TH1D>(new TH1D{("zResidualPull"), "Pull of z-residual;z-residual pull; Hits", 100,-4,4});

	pixelHitMap=unique_ptr<TH2D>(new TH2D{"pixelHitMap", "Hitmap by pixel;Columns;Rows", 256,0, 256, 256,0,256});

	xRotation=unique_ptr<TH1D>(new TH1D{"xRotation", "x-rotation;x-rotation [rad.]; Hits", 100,-0.2,0.2});
	yRotation=unique_ptr<TH1D>(new TH1D{"yRotation", "y-rotation;y-rotation [rad.]; Hits", 100,-0.2,0.2});
	zRotation=unique_ptr<TH1D>(new TH1D{"zRotation", "z-rotation;z-rotation [rad.]; Hits", 100,-0.2,0.2});

	xzShear=unique_ptr<TH1D>(new TH1D{"xzShear", "xz-shear;xz-shear; Hits", 100,-0.5,0.5});


	zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	xResidualByToT=std::unique_ptr<TH2D>(new TH2D{"xResidualByToT", "x-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	zResidualByDriftTime=std::unique_ptr<TH2D>(new TH2D{"zResidualByDriftTime", "z-residuals as a function of drift time;Drift time [#mus];z-residual [mm]", int(0.8/ToABinWidth),-0.39999,0.4,50,-2,2});

	xByzForShear=std::unique_ptr<TH2D>(new TH2D{"xByzForShear", "x-residual by z;z-position [mm]; x-residual [mm]", 100,-2,12, 100,-1,1});

	auto rangeGlobal=rangeFromCorners( align.getAllChipCorners() );
//	auto rangeQuad=rangeFromCorners( align.getAllChipCornersQuad() ); //automatic range
	PositionRange rangeQuad{0,29.04,-14.55,25.05};

	global.createHistograms(rangeGlobal);
	local.createHistograms(rangeQuad);
	locExp.createHistograms(rangeQuad);
	xResidualByPixel=std::unique_ptr<TProfile2D>(new TProfile2D{"xResidualByPixel", "mean x-residual by pixel;Columns;Rows", 256,0,256,256,0,256});
	zResidualByPixel=std::unique_ptr<TProfile2D>(new TProfile2D{"zResidualByPixel", "mean z-residual by pixel;Columns;Rows", 256,0,256,256,0,256});

	xResidualByPixel->SetMinimum(-0.1), xResidualByPixel->SetMaximum(0.1);
	zResidualByPixel->SetMinimum(-0.1), zResidualByPixel->SetMaximum(0.1);

	startDir->cd();
}

void ChipHistogrammer::residualHistograms::createHistograms(const PositionRange& range) {
	using namespace std;
	const int nBins=160, nBinsZ=200;

	auto startDir=gDirectory;
	if(!startDir) std::cout<<"current directory is not valid!\n";
	auto newDir=startDir->mkdir( name.substr(1,std::string::npos).c_str() );
	if(!newDir) std::cout<<"could not make directory: " <<name<<"\n";
	newDir->cd();

	xResidual=unique_ptr<TH1D>(new TH1D{("xResidual"), "x-residual;x-residual [mm]; Hits", nBins,-2,2});
	yResidual=unique_ptr<TH1D>(new TH1D{("yResidual"), "y-residual;y-residual [mm]; Hits", nBins,-2,2});
	zResidual=unique_ptr<TH1D>(new TH1D{("zResidual"), "z-residual;z-residual [mm]; Hits", nBins,-4,4});

	double zmin=-10, zmax=50;
	xResidualByz=std::unique_ptr<TH2D>(new TH2D{("xResidualByz"), "x-residuals as a function of drift distance;Drift distance [mm];x-residual [mm]", nBinsZ,zmin,zmax,nBins,-2,2});
	yResidualByz=std::unique_ptr<TH2D>(new TH2D{("yResidualByz"), "y-residuals as a function of drift distance;Drift distance [mm];y-residual [mm]", nBinsZ,zmin,zmax,nBins,-2,2});
	zResidualByz=std::unique_ptr<TH2D>(new TH2D{("zResidualByz"), "z-residuals as a function of drift distance;Drift distance [mm];z-residual [mm]", nBinsZ,zmin,zmax,nBins,-2,2});
	zResidualByzByToT=std::unique_ptr<TH3D>(new TH3D{("zResidualByzByToT"), "z-residuals as a function of drift distance;Drift distance [mm];z-residual [mm];ToT [#mus]", nBinsZ,-1,14,nBins,-2,2,100,0,2.5});

	zHit = unique_ptr<TH1D>(new TH1D{ ("zHit"), "z-position of hits;z [mm]; Hits", 200,zmin,zmax });
	double xmax=range.xmax, xmin=range.xmin, ymin=range.ymin, ymax=range.ymax;
	positionHitMap=unique_ptr<TH2D>(new TH2D{ ("positionHitMap"), "Hitmap with positions;x-position [mm];y-position [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});
	XZHitMap=unique_ptr<TH2D>(new TH2D{ ("XZHitMap"), "Hitmap;x-position [mm];z-position [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((zmax-zmin)/0.2),zmin, zmax});
	YZHitMap=unique_ptr<TH2D>(new TH2D{ ("YZHitMap"), "Hitmap;y-position [mm];z-position [mm]",
		int((ymax-ymin)/.055),ymin,ymax,int((zmax-zmin)/0.2),zmin, zmax});


	xResidualByPosition=std::unique_ptr<TProfile2D>(new TProfile2D{ ("xResidualByPosition"), "mean x-residual by position;x-position [mm]; y-position [mm]; mean x-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});
	zResidualByPosition=std::unique_ptr<TProfile2D>(new TProfile2D{ ("zResidualByPosition"), "mean z-residual by position;x-position [mm]; y-position [mm]; mean z-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});
	for(auto prof2d : {&xResidualByPosition, &zResidualByPosition} ) {
		(*prof2d)->SetMinimum(-0.1), (*prof2d)->SetMaximum(0.1);
	}

	xResidualByXZ=std::unique_ptr<TProfile2D>(new TProfile2D{ ("xResidualByXZ"), "mean x-residual by XZ;x-position [mm]; z-position [mm]; mean x-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,200,zmin, zmax});
	zResidualByXZ=std::unique_ptr<TProfile2D>(new TProfile2D{ ("zResidualByXZ"), "mean z-residual by XZ;x-position [mm]; z-position [mm]; mean z-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,200,zmin, zmax});
	for(auto prof2d : {&xResidualByXZ, &zResidualByXZ} ) {
		(*prof2d)->SetMinimum(-0.2), (*prof2d)->SetMaximum(0.2);
	}

	startDir->cd();
}

void ChipHistogrammer::fillHit(const PositionHit& h) {
	driftTime->Fill(h.driftTime/4096.*25E-3);
	ToT->Fill(h.ToT*25E-3);
	pixelHitMap->Fill(h.column, h.row);
	zResidualByToTCorrected->Fill(h.ToT*25E-3, h.residual.z);
	xResidualByToT->Fill(h.ToT*25E-3, h.residual.x);
	zResidualByDriftTime->Fill(h.driftTime/4096.*25E-3,h.residual.z);

	xResidualPull->Fill(h.residual.x/h.error.x);
	yResidualPull->Fill(h.residual.y/h.error.y);
	zResidualPull->Fill(h.residual.z/h.error.z);

	global.fill(h.position, h.residual, h.ToT);
	xResidualByPixel->Fill(h.column,h.row, h.residual.x);
	zResidualByPixel->Fill(h.column,h.row, h.residual.z);

	hitsCounter++;
}

void ChipHistogrammer::fillRotation(const TVector3& position, const TVector3& residual, const TVector3& COM) { //h in timepix coordinates, COM
	//todo: do proper 3d rotation alignment (not just small angle approx)
	auto d=position-COM;

	//x
	double phix=residual.z()/d.y();
	xRotation->Fill(phix, d.y()*d.y());

	//y
	double perp2y=d.Perp2({0,1,0});
	double phiy=(d.x()*residual.z()-d.z()*residual.x())/perp2y;
	yRotation->Fill(phiy, perp2y);
	//for the moment, until error are included, use only x residuals in y rotation!
//	double phiy=-h.residual.x/d.z();
//	yRotation->Fill(phiy,d.z()*d.z());

	xzShear->Fill( residual.x()/d.z(), d.z()*d.z() );
	xByzForShear->Fill(d.z(), residual.x());

	//z
//	double phi=(d.y()*h.residual.x-d.x()*h.residual.y)/d.Perp2();
//	double weight=d.Perp2();
	double phi=residual.x()/d.y();
	double weight=d.y()*d.y();
	zRotation->Fill(phi,weight);
}
void ChipHistogrammer::fillEvent() {
	if(!hitsCounter) return;
	nHits->Fill(hitsCounter);
	hitsCounter=0;
}


QuadTrackFitter::QuadTrackFitter(std::string fileName) :
		reader("data", fileName) {}

QuadTrackFitter::~QuadTrackFitter() {
}

void QuadTrackFitter::Loop(std::string outputFile,const Alignment& alignment) {

	TFile file(outputFile.c_str(), "RECREATE");

	Vec3 averageHitPosition;
	std::vector<int> nHits(4);
	int nHitsPassedTotal{0};

	TTree fitResults("fitResults", "fitResults");
	fitResults.Branch("hits", "std::vector<PositionHit>", &posHits);
	fitResults.Branch("laser", "Vec3", &laser);
	fitResults.Branch("hitAverage", "Vec3", &averageHitPosition);
	fitResults.Branch("nHitsPassedChip", "std::vector<int>", &nHits);
	fitResults.Branch("nHitsPassed", &nHitsPassedTotal);

  auto zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});
  auto zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});

	std::array<ChipHistogrammer, nChips> hists{{
		{"chip1", alignment},
		{"chip2", alignment},
		{"chip3", alignment},
		{"chip4", alignment} }};
	ChipHistogrammer quad{"quad", alignment};
	reader.tree->GetEntry();

	 do {
		const auto driftSpeed=alignment.driftSpeed.value/4096*25; //in units of mm/ticks
		posHits=convertHitsQuad(reader.chip, driftSpeed);
		laser=Vec3(reader.laser->x, reader.laser->y, reader.laser->z);//mirror z-axis, also mirror x-axis to keep cooridnate system right handed

		if(posHits.size()>4000) continue;//this is probably a discharge event
		if(not selectLaserPoint(laser)) continue;

		nHits=std::vector<int>(4);
		for(auto& h : posHits) {
			//apply alignment
			h=alignment.transform(h);
			h=alignment.timeWalk.correct(h);
//			if(h.position.z<alignment.hitErrors.z0) h.flag=PositionHit::Flag::smallz;
			h.error=alignment.hitErrors.hitError(h.position.z);

			h.calculateResidual(laser);
			h=flagResidual(h, {2,4,2});
//			h=flagResidualPull(h, {2,2,3});

			if(h.flag==PositionHit::Flag::valid) nHits[h.chip]++;
		}
		nHitsPassedTotal=std::accumulate(nHits.begin(), nHits.end(),0);

		averageHitPosition={0,0,0};
		for(auto& h : posHits) {

			if(h.flag==PositionHit::Flag::highResidualxy) continue;
			hists[h.chip].fillTimewalk(h.residual.z, h.ToT, alignment.timeWalk);
			quad.fillTimewalk(h.residual.z, h.ToT, alignment.timeWalk);
			zResidualByToTCorrected->Fill(h.ToT*25E-3,h.residual.z);
			zResidualByToT->Fill(h.ToT*25E-3,h.residual.z+alignment.timeWalk.getCorrection(h.ToT));

			if(h.flag!=PositionHit::Flag::valid) continue;
			hists[h.chip].fillHit(h);

			quad.fillHit(h);
			quad.fillRotation(h, alignment.quad.getShiftedCOM());
			auto localPosition=alignment.quad.rotateAndShiftBack(h.position);
			auto localResidual=alignment.quad.rotateBack(h.residual, {0,0,0} );
			auto expectedPosition = h.position-h.residual;
			auto localExpectedPosition=expectedPosition;
			localExpectedPosition=alignment.quad.rotateAndShiftBack(localExpectedPosition);

			auto chipFrameResidual=alignment.chips[h.chip].rotateBack(localResidual, {0,0,0});
			if(h.chip==0 || h.chip==3) {
				chipFrameResidual.x*=-1; chipFrameResidual.y*=-1;
			}
			auto chipFrameExpectedPosition=localExpectedPosition;
			chipFrameExpectedPosition=alignment.chips[h.chip].rotateAndShiftBack(chipFrameExpectedPosition);
			hists[h.chip].fillRotation(chipFrameExpectedPosition, chipFrameResidual, alignment.chips[h.chip].COM);

			hists[h.chip].local.fill( localPosition, localResidual, h.ToT );
			quad.local.fill( localPosition, localResidual, h.ToT );

//			hists[h.chip].locExp.fill( localExpectedPosition, localResidual, h.ToT );
			hists[h.chip].locExp.fill( chipFrameExpectedPosition, chipFrameResidual, h.ToT );
			quad.locExp.fill( localExpectedPosition, localResidual, h.ToT );

			averageHitPosition=averageHitPosition+h.position;
		}
		if(nHitsPassedTotal) averageHitPosition=1./nHitsPassedTotal*averageHitPosition;

		for(auto& h : hists) {
			h.fillEvent();
		}
		quad.fillEvent();

		fitResults.Fill();

		const bool draw=false;
		if(draw && laser.x>=4 and not (reader.entryNumber%200)) {
//			drawCluster2DPixel(posHits);

			static auto setupForDrawing=simpleDetectorFromChipCorners(alignment.getAllChipCorners());
			setupForDrawing.minz=std::min(setupForDrawing.minz,10.0);
			setupForDrawing.maxz=std::max(setupForDrawing.maxz,50.0);;

			//2D
			drawCluster2D(posHits,setupForDrawing);
			alignment.drawChipEdges();
			drawLaserTrack2D(laser.x, setupForDrawing.ymin(), setupForDrawing.ymax());
			gPad->Update();

			//3D
			drawCluster(posHits,setupForDrawing);
			drawQuadOutline(alignment, setupForDrawing.zmax() );
			drawLaserTrack3D(laser.x, laser.z, setupForDrawing.ymin(), setupForDrawing.ymax());
			gPad->GetPrimitive("eventDisplayLegend")->Draw(); //move legend to top again
			gPad->Update();

			if( processDrawSignals() ) return;
		}

	} while(reader.getNextEntry());

	file.Write();
}
