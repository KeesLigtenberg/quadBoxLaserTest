/*
 * drawFunctions.h
 *
 *  Created on: Jun 12, 2019
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_DRAWFUNCTIONS_H_
#define LASERDATAFITTER_DRAWFUNCTIONS_H_

#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TView.h"
#include "TPaletteAxis.h"
#include "TPaveLabel.h"
#include "TPolyLine3D.h"
#include "TLegend.h"

template<class T > //=std::list<HitCluster>
inline void drawCluster2D(const T& cluster, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
//	pointTree.Fill(); //fill one with zero ToT to set scale

	for(auto& iHit : cluster ) {
		h=iHit;
//		if(h.flag==PositionHit::Flag::valid)
			pointTree.Fill();
	}
	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(7);
	gStyle->SetOptTitle(0);
	double totAxis=1.6;
	static TCanvas* canv=new TCanvas("cluster_display2D", "Event display 2D", 600,800);
	canv->cd();
	pointTree.Draw( "h.position.y:h.position.x" , "", ""); //ToT to microseconds

	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("y-axis (beam direction) [mm]");
	yaxis->SetLimits(detector.ymin(),detector.ymax());
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetLimits(detector.xmin(),detector.xmax());
	xaxis->SetTitle("x-axis [mm]");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.1);
	for(auto axis : {xaxis, yaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
//	axisObject->SetMaximum(totAxis);
//	axisObject->SetMinimum(0);

	axisObject->Draw();
	pointTree.Draw( "h.position.y:h.position.x" , "", "same");
	gPad->SetMargin(0.15,0.1,0.1,0.05);//l r b t
	gPad->Update();

//
//	TLegend* legend= new TLegend( 0.6, 0.8, 0.95,0.95 );
//	legend->SetName("eventDisplayLegend");
//	legend->AddEntry(axisObject, "Timepix hits", "p");
//	axisObject->SetLineColor(kOrange+7);
//	axisObject->SetLineWidth(2);
//	legend->AddEntry(axisObject, "Telescope track", "l");
//	legend->Draw();

	gPad->Update();

}

template<class T > //=std::list<HitCluster>
inline void drawCluster2DPixel(const T& cluster) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
//	pointTree.Fill(); //fill one with zero ToT to set scale

	TProfile2D prof("pixelProfile", "pixelProfile", 512,0,512,512,0,512);

	for(auto& iHit : cluster ) {
		h=iHit;
//		std::cout<<int(h.ToT)<<"\n";
//		if(h.flag==PositionHit::Flag::valid)
//			pointTree.Fill();
				prof.Fill(
						(h.chip>=2)*256+(h.chip==0||h.chip==3)*(256-2*h.column)+h.column,
						(h.chip==1||h.chip==2)+255+(h.chip==0||h.chip==3)*-2*h.row+h.row,
						h.driftTime/4096.*25E-3);//choose between driftTime,
//						h.ToT*0.025); //and ToT
	}
//	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(7);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kRainBow);
	double totAxis=5;
	static TCanvas* canv=new TCanvas("cluster_display2DPixel", "Event display 2D per pixel", 900,800);
	canv->cd();
//	pointTree.Draw( "h.ToT"
//			":(h.chip==1||h.chip==2)+255+(h.chip==0||h.chip==3)*-2*h.row+h.row"
//			":(h.chip>=2)*256+(h.chip==0||h.chip==3)*(256-2*h.column)+h.column"
//			">>hpixel(512,0,512,512,0,512)" , "", "profcolz");

	auto axisObject=&prof; //dynamic_cast<TH1*>( gPad->GetPrimitive("hpixel") );
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("Rows");
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetTitle("Columns");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.1);
	for(auto axis : {xaxis, yaxis} ) {
		axis->SetTitleSize(0.045);
		axis->SetLabelSize(0.045);
	}
	axisObject->SetMaximum(totAxis);
	axisObject->SetMinimum(0);

	axisObject->DrawCopy("colz0");
	//pointTree.Draw( "h.position.y:h.position.x" , "", "same");
	gPad->SetMargin(0.15,0.15,0.15,0.05);//l r b t
	gPad->Update();

	gPad->Update();

}


//todo: move all draw functions away from houghtransformer
inline bool processDrawSignals() {
	static bool printedInfo=false;
	static bool pdfOpen=false;
	if(not printedInfo) { std::cout<<"<return> to continue, 'q' to break\n"; printedInfo=true; };

	char userSignal = std::cin.get();
//	std::cout<<"signal = "<<userSignal<<"\n";
	if (userSignal == 'q') {
		if(pdfOpen) {gPad->Print("eventDisplays.pdf]"); pdfOpen=false; }//close pdf
		return true;
	} else if (userSignal == 'l') {
		while (!gSystem->ProcessEvents()) {
			gSystem->Sleep(50);
		}
		return true;
	} else if (userSignal == 'w') {
		//rotate and write as animated gif!
		double phiView = 55;
		for (int thetaView = 0; thetaView < 360; thetaView += 2) {
			gPad->GetView()->RotateView(thetaView, phiView);
			gPad->Modified();
			gPad->Update();
			gPad->Print( thetaView == 358 ?	"eventAnimation.gif++5++" : "eventAnimation.gif+5");
//			gPad->Print(( "eventAnimation"+to_string(thetaView/2)+".png" ).c_str());
		}
	} else if (userSignal == 'a') {
		if(not pdfOpen) {
			gPad->Print("eventDisplays.pdf(");
			pdfOpen=true;
		} else {
			gPad->Print("eventDisplays.pdf");
		}
	}
	return false;
}

void drawLaserTrack2D(double xlaser, double zmin, double zmax, Color_t color=kOrange+7) { //global coordinates!
	const int npoints=2;
	double x[npoints] = {xlaser, xlaser};
	double z[npoints] = { zmin, zmax};
	TPolyLine l( npoints, x, z );
	l.SetLineColor(color);
	l.SetLineWidth(2);
	l.DrawClone();
}
void drawLaserTrack3D(double xlaser, double ylaser, double zmin, double zmax, Color_t color=kOrange+7) {  //global coordinates!
	const int npoints=2;
	double x[npoints] = { xlaser, xlaser};
	double y[npoints] = { ylaser, ylaser};
	double z[npoints] = { zmin, zmax};
	TPolyLine3D l( npoints, x, z, y );
	l.SetLineColor(color);
	l.SetLineWidth(2);
	l.DrawClone();
}

template<class T>
void drawCluster(const T& cluster, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
	pointTree.Fill(); //fill one with zero ToT to set scale

	for(auto& iHit : cluster ) {
		h=iHit;
//		std::cout<<int(h.ToT)<<"\n";
//		if(h.flag==PositionHit::Flag::shiftedTrigger) continue;
//		if(h.flag==PositionHit::Flag::valid)
			pointTree.Fill();
	}
	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(20);
	gStyle->SetOptTitle(0);
	double totAxis=1.6;
	static TCanvas* canv=new TCanvas("cluster_display", "Event display", 800,600);
	canv->cd();
//	pointTree.Draw( std::string("h.position.z:h.position.y:h.position.x:h.nShiftedTrigger").c_str() , "", "*colz"); //ToT to microseconds
	pointTree.Draw( ("h.position.z:h.position.y:h.position.x:TMath::Min(h.ToT*0.025, "+std::to_string(totAxis)+")").c_str() , "", "*colz"); //ToT to microseconds
//	pointTree.Draw( std::string("h.position.z:h.residual.y:h.residual.x:h.chip").c_str() , "", "*colz"); //ToT to microseconds

	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	if(!axisObject) {std::cout<<"could not get axis object!?\n"; throw 1;}
	auto zaxis=axisObject->GetZaxis();
	if(!zaxis) {std::cout<<"could not get zaxis!?\n"; throw 1;}
	zaxis->SetLimits(detector.zmin(),detector.zmax());
	zaxis->SetTitle("z-axis (drift direction) [mm]") ;
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("y-axis (beam direction) [mm]");
	yaxis->SetLimits(detector.ymin(),detector.ymax());
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetLimits(detector.xmin(),detector.xmax());
	xaxis->SetTitle("x-axis [mm]");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.3);
	for(auto axis : {xaxis, yaxis,zaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
//	axisObject->SetMaximum(totAxis);
//	axisObject->SetMinimum(0);

	axisObject->Draw("colz");
	gPad->SetMargin(0.1,0.175,0.15,0.1);//l r b t
	gPad->Update();

	TPaletteAxis* palette= dynamic_cast<TPaletteAxis*>(gPad->GetPrimitive("palette"));
	if(!palette) throw "could not find paletteAxis!";
	palette->SetX1NDC(0.85);
	palette->SetX2NDC(0.90);
	palette->SetY2NDC(0.74);
	palette->SetY1NDC(0.1);
	//draw TPaveText over Palette axis title
	auto paletteAxisLabel = new TPaveLabel(0.96,0.1,1,0.75, "ToT [#mus]", "NDC");
	paletteAxisLabel->SetFillColor(kWhite);
	paletteAxisLabel->SetBorderSize(0);
	paletteAxisLabel->SetTextAngle(90);
	paletteAxisLabel->SetTextSize(0.08);
	paletteAxisLabel->SetTextFont(42);
	paletteAxisLabel->SetTextAlign(kHAlignCenter+kVAlignCenter);
	paletteAxisLabel->Draw();

	double theta=-20 /*-20*/,phi=60 /*10*/;
	//	std::cout<<"give angles!"<<std::endl;
	//	std::cin>>theta>>phi;
	gPad->GetView()->RotateView(theta, phi);

	TLegend* legend= new TLegend( 0.6, 0.8, 0.95,0.95 );
	legend->SetName("eventDisplayLegend");
	legend->AddEntry(axisObject, "Timepix hits", "p");
	axisObject->SetLineColor(kOrange+7);
	axisObject->SetLineWidth(2);
	legend->AddEntry(axisObject, "Telescope track", "l");
	legend->Draw();

	gPad->Update();

}

void drawChipFromCorners(std::array<TVector3, 4> corners) {
	const int npoints=5;
	double x[npoints] = {};
	double y[npoints] = {};
	double z[npoints] = {};
	for(int i=0; i<npoints; i++) {
		x[i]=corners[i%4].x();
		y[i]=corners[i%4].y();
		z[i]=corners[i%4].z();
	}
	TPolyLine3D l( npoints, x, y, z );
	l.SetLineColor(kAzure-8);
	l.SetLineWidth(2);
	l.DrawClone();
}

void drawQuadOutline(const Alignment& align, double z) {
	for(unsigned i=0; i<align.chips.size(); i++) {
		auto corners=align.getChipCorners(i);
		for(auto& c: corners) c.SetZ( z );
		drawChipFromCorners(corners);
	}
}

#endif /* LASERDATAFITTER_DRAWFUNCTIONS_H_ */
