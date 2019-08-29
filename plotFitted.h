//root macro
#include <cstdlib>
#include <vector>

#include "TPolyLine.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH3.h"
#include "TF2.h"
#include "TText.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "/user/cligtenb/rootmacros/AllCombiner.h"
#include "/user/cligtenb/rootmacros/histogramOperations.h"
#include "/user/cligtenb/rootmacros/StatsWrapper.h"
#include "/user/cligtenb/rootmacros/CombineHistogramsFromTree.h"
#include "/user/cligtenb/rootmacros/CombineHistogramsFromFiles.h"
#include "/user/cligtenb/rootmacros/getMeanFromFit.h"

//#include "../rootmacros/getObjectFromFile.h"
//#include "../rootmacros/getHistFromTree.h"
//#include "../rootmacros/AllCombiner.h"
//#include "../rootmacros/histogramOperations.h"
//#include "../rootmacros/StatsWrapper.h"


#include "laserDataFitter/Alignment.h"
//#include "TrackCombiner/deformationCorrection.h"

using namespace std;

const int nChips=4;
std::vector<std::string> chipDirectories={"chip1", "chip2","chip3","chip4"};
std::vector<int> chipMapping={2,1,3,0};
std::vector<int> chipReverseMapping={4,2,1,3};//for canv->cd

void drawChipEdges(std::string alignFile="align.dat") {
	Alignment	alignment(alignFile);
	alignment.drawChipEdges();
	gStyle->SetOptStat(0);
}

void drawChipEdgesLocal(std::string alignFile="align.dat") {
	Alignment	alignment(alignFile);
	alignment.drawChipEdges(false);
	gStyle->SetOptStat(0);
}

void combineHistogramsForChips(std::string histName="nHits", std::string fileName="combinedFit.root", bool addQuad=false) {
	HistogramCombiner combination(histName+"combination");

	if(addQuad) {
		combination.setStyle(3);
		auto quadHist=getObjectFromFile<TH1>("quad/"+histName, fileName);
		if(histName=="nHits") {
			quadHist->Rebin(2);
			quadHist=convertToTH1D(quadHist);
		}
		combination.add( quadHist , "All chips" );
	}

	for(int i=0; i<4; i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
		if(histName=="nHits") {
			hist->Rebin(2);
			hist=convertToTH1D(hist);
		}
//		hist->Draw(); gPad->Update(); cin.get();
		combination.add( hist , "chip "+std::to_string(i+1) );
	}
	combination.normalise();
	combination.createCombined();
}
void combineDrawHistogramsForChips(std::string expression, std::string cut="1", std::string drawOption="", std::string histName="Hist", std::string fileName="combinedFit.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	HistogramCombiner combination(histName+"Combination");
	for(int i=0;i<4;i++) {
		auto hist=getHistFromTree(tree, expression, cut+" && chip=="+to_string(i), "chip"+to_string(i)+histName, drawOption.c_str(),1E6);
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
		combination.add(hist);
	}
	combination.createCombined();
}

void combineHistogramsForChipsOnPads(std::string histName="nHits", std::string fileName="fitted.root") {
	auto canv=new TCanvas(("canv"+histName).c_str(), "combined canvas", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
		canv->cd(i+1);
		hist->Draw("colz");
		getMeanFromGausFitAroundMean(*hist);
	}
}

//draw four histograms on 4 pads in one canvas
void combineDrawPadsForChips(std::string expression, std::string cut, std::string drawOption, std::string histName="Hist", std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto hist=getHistFromTree(tree, expression, cut+" && chip=="+to_string(i), "chip"+to_string(i)+histName, drawOption.c_str());
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
	}
}

//only accapted hits
void plotHitMap(std::string file="combinedFit.root", std::string alignFile="align.dat", std::string hitMapName="global/positionHitMap") {
	new TCanvas("hitmap", "hitmap", 850,1000)	;

	//hits
	auto first=true;

//	for(const auto& dir:chipDirectories)
	std::string dir="quad";
	{
		auto hist=getObjectFromFile<TH2>(dir+"/"+hitMapName, file);
		if(first) {
			hist->SetTitle(";x-position [mm];y-position [mm]; Hits");
			hist->Draw("colz");
			gPad->SetMargin(0.1,0.15,0.1,0.05);
			gPad->SetTicks(1,1);

			hist->GetYaxis()->SetTitleOffset(1.3);
			hist->GetXaxis()->SetTitleOffset(1.1);
			hist->GetZaxis()->SetTitleOffset(1.3);
			
		} else {
			hist->Draw("colsame");
		}
		hist->Draw( first ? (first=false, "colz") : "colsame" );
	}

	//chip edges
	if(getSubStringBefore(hitMapName, "/")=="local")
		drawChipEdgesLocal(alignFile);
	else
		drawChipEdges(alignFile);
}

//this is all hits
void plotHitMapFromDraw(std::string file="fitted.root", std::string alignFile="align.dat") {

	new TCanvas("hitmap", "hitmap", 850,1000)	;
	auto tree=getObjectFromFile<TTree>("fitResults",file);
	auto hist=getHistFromTree(*tree, "position.y:position.x", "", "h(564,-5,27,745,173,214)", "colz",1E6);

	hist->SetTitle(";x-position [mm];y-position [mm]; Hits");
	gPad->SetMargin(0.1,0.15,0.1,0.05);
	gPad->SetTicks(1,1);

	hist->GetYaxis()->SetTitleOffset(1.3);
	hist->GetXaxis()->SetTitleOffset(1.1);
	hist->GetZaxis()->SetTitleOffset(1.3);

	drawChipEdges(alignFile);

}

void combineHistogramsWithMean(std::string histName="xResidual", std::string fileName="combinedFit.root") {
	HistogramCombiner combination(histName+"combination");
	for(int i=0; i<4; i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		auto meanstr=std::to_string(int(1000*hist->GetMean()));
		combination.add( hist , "chip "+std::to_string(i+1)+" (mean "+meanstr+" #mum)" );
	}
	//combination.normalise();
	combination.createCombined();
}

void drawToT(std::string fileName="combinedFit.root") {
	combineDrawHistogramsForChips("ToT*25E-3", "fabs(ToT*25E-3)<2.5", "", "hist(100,0,2.5)", fileName);
}

void combineColumnToACorrection(std::string histName="zResidualByPixel", std::string fileName="combinedFit.root") {
	HistogramCombiner combination(histName+"combination");
	for(int i=0; i<4; i++) {
		auto prof=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		auto hist=prof->ProfileX();
		combination.add( hist , "chip "+std::to_string(i+1) );
	}
	//combination.normalise();
	combination.createCombined();

	new TCanvas("quad","quad", 500,500);
	auto prof=getObjectFromFile<TProfile2D>("quad/"+histName, fileName);
	prof->ProfileX()->Draw();
}

void plotDeformations(std::string x="x", std::string file="fitted.root", std::string alignFile="../align.dat") {
	new TCanvas(("deformations"+x).c_str(), ("deformations"+x).c_str(), 900,1000)	;
	
	auto tree=getObjectFromFile<TTree>("fitResults", file);
	//auto deformation=getHistFromTree(tree, "hitAverage."+x+"-laser."+x+":laser.y:laser.x", "hitAverage.x>0", "deformation(29,9.5,38.5,35, 1.5, 40.5)", "profcolz");
	auto deformation=getHistFromTree(tree, "residual."+x+":laser.y:laser.x", "flag>0 && hitAverage.x>0", "deformation(29,11.5,40.5,35, 1.5, 40.5)", "profcolz");
	deformation->SetBins(29,11.5,40.5,39, 1.5, 40.5);
	setMinMax((TProfile2D*) deformation, -0.2,0.2);
	deformation->SetTitle((";x [mm];y [mm];Mean "+x+"-residual [mm]").c_str());
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	deformation->GetYaxis()->SetTitleOffset(1.3);
	deformation->GetXaxis()->SetTitleOffset(1.1);
	deformation->GetZaxis()->SetTitleOffset(1.6);

	drawChipEdges(alignFile);
}
void plotDeformationsXYZ(std::string file="fitted.root", std::string alignFile="../align.dat") {
	plotDeformations("x", file,alignFile);
	plotDeformations("y", file,alignFile);
	plotDeformations("z", file,alignFile);
}
//draw deformations as pixel coordinates
void plotDeformationsPixels(std::string x="x", std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas(("canvasDeformationPixels"+x).c_str(), "combined canvas", 1200,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto chipHist=new TProfile2D(("chipHist"+to_string(i)).c_str(), "-residuals per pixel", 14,0,256, 14,0,256);
		getHistFromTree(tree, "residual."+x+":row:column", "flag>0 && hitAverage.x>0 && chip=="+to_string(i), "chipHist"+to_string(i), "profcolz0");
		//chipHist->SetBins(256,0,256, 256,0,256);
		chipHist->SetMinimum(-0.2);
		chipHist->SetMaximum(0.2);
		//setMinMax((TProfile2D*) chipHist,-0.2,0.2);
		chipHist->SetTitle(("Chip "+to_string(i+1)+";Columns;Rows;Mean "+x+"-residual [mm]").c_str());
		gPad->SetMargin(0.1,0.2,0.1,0.1);
		gPad->SetTicks(1,1);
		chipHist->GetYaxis()->SetTitleOffset(1.3);
		chipHist->GetXaxis()->SetTitleOffset(1.1);
		chipHist->GetZaxis()->SetTitleOffset(1.6);
		
		//chipHist->GetXaxis()->SetLabelSize(0.05);
		chipHist->GetXaxis()->SetNdivisions(408, false);
		//chipHist->GetYaxis()->SetLabelSize(0.05);
		chipHist->GetYaxis()->SetNdivisions(408, false);
	}
	gStyle->SetOptStat(0);
}

//draw deformations as pixel coordinates
void plotToTPerPixel(std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas("canvas", "combined canvas", 1200,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto chipHist=new TProfile2D(("chipHist"+to_string(i)).c_str(), "-residuals per pixel", 32,0,256, 32,0,256);
		getHistFromTree(tree, "ToT/40.:row:column", "flag>0 && hitAverage.x>0 && chip=="+to_string(i), "chipHist"+to_string(i), "profcolz0");
		//chipHist->SetBins(256,0,256, 256,0,256);
		chipHist->SetMinimum(0.3);
		chipHist->SetMaximum(0.7);
		//setMinMax((TProfile2D*) chipHist,0,2);
		chipHist->SetTitle(("Chip "+to_string(i+1)+";Columns;Rows;Mean ToT [#mus]").c_str());
		gPad->SetMargin(0.1,0.2,0.1,0.1);
		gPad->SetTicks(1,1);
		chipHist->GetYaxis()->SetTitleOffset(1.3);
		chipHist->GetXaxis()->SetTitleOffset(1.1);
		chipHist->GetZaxis()->SetTitleOffset(1.6);
		
		//chipHist->GetXaxis()->SetLabelSize(0.05);
		chipHist->GetXaxis()->SetNdivisions(408, false);
		//chipHist->GetYaxis()->SetLabelSize(0.05);
		chipHist->GetYaxis()->SetNdivisions(408, false);
	}
	gStyle->SetOptStat(0);
}

void plotHitAverage(std::string file="fitted.root") {
	new TCanvas("hitAverage", "hitAverage", 750,1000)	;
	TTree* tree=getObjectFromFile<TTree>("fitResults",file);

	//hit positions
	TGraph* hitGraph=getGraphFromTree(*tree,"hitAverage.y:hitAverage.x", "hitAverage.x>0 && nHitsPassed>3");
	//hitGraph->SetMarkerStyle(7);
	hitGraph->SetTitle(";x [mm];y [mm]");
	hitGraph->Draw("AP");

	//chip edges
	drawChipEdges();
}

void plotResidualsTimeWalk( std::string filename="fitted.root" ) {

	auto uncorrected=getObjectFromFile<TH2>("quad/zResidualByToT", filename);
	TH1* uncorredtedHist=uncorrected->ProjectionY();
	auto corrected=getObjectFromFile<TH2>("quad/zResidualByToTCorrected", filename);
	TH1* correctedHist=corrected->ProjectionY();
	uncorredtedHist=makeXShifted(uncorredtedHist, -uncorredtedHist->GetMean() );

	HistogramCombiner combination("zResiduals time walk");
	combination.add(uncorredtedHist, "z-residuals without time walk correction;z-residual [mm]; Hits");
	combination.add(correctedHist, "z-residuals with time walk correction");
//	combination.normalise();
	combination.setNcolumn(1);
	combination.createCombined();
	TGaxis::SetMaxDigits(4);
}

TH1* plotFittedTimeWalk(TH2* uncorrected) {
	uncorrected->FitSlicesY();
	auto means=dynamic_cast<TH1*>( gDirectory->Get("zResidualByToT_1") );
	if(!means) { cerr<<"failed to retrieve result from fitslicesy()\n"; return nullptr; };

	//rebin means
	if (false) {
		auto hist=means;
		std::vector<double> lowerEdges, contents, errors;
		lowerEdges.push_back( hist->GetXaxis()->GetBinLowEdge(1));
		for(int bin=1; bin<=hist->GetNbinsX();) {			
			double lowEdge=hist->GetXaxis()->GetBinLowEdge(bin);
			int nBinsPerBin=  lowEdge<1.5 ? 1 : lowEdge<2 ? 2 : 4;

			double content=0, error2=0;
			for(int i=0; i<nBinsPerBin; i++) {
				content+=hist->GetBinContent(bin);
				error2+=hist->GetBinError(bin)*hist->GetBinError(bin);
				bin++;
			}
			contents.push_back(content/nBinsPerBin);
			errors.push_back(sqrt(error2)/nBinsPerBin);
			lowerEdges.push_back(hist->GetXaxis()->GetBinLowEdge(bin));

		}
		auto rebinned=new TH1D("rebinned", "rebinned", lowerEdges.size()-1, lowerEdges.data() );
		for(int i=0; i<contents.size(); i++) {
			rebinned->SetBinContent(i+1,contents[i]);
			rebinned->SetBinError(i+1,errors[i]);
		}
		means=rebinned;
	}
	

	double minToT=0.10;
	auto simple=new TF1("2 parameters", "[c_{1}]/(x+[t_{0}])+[offset]", minToT,2.5);
	simple->SetLineWidth(1);
	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);
	double offset=simple->GetParameter("offset");
	means=shiftY(means, -simple->GetParameter("offset") );

	new TCanvas("fittedTW", "fitted TW", 800,600);
	means->Draw();	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);


	means->GetYaxis()->SetTitle("Mean z-residual [mm]");
	means->GetYaxis()->SetTitleOffset(1.3);
	means->GetYaxis()->SetTitleSize(0.05);
	means->GetYaxis()->SetLabelSize(0.05);
	means->GetYaxis()->SetRangeUser(0,3);
	means->GetXaxis()->SetTitleOffset(1.3);
	means->GetXaxis()->SetTitleSize(0.05);
	means->GetXaxis()->SetLabelSize(0.05);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	simple->SetLineWidth(1);
	simple->SetNpx(1000);
	gPad->SetTicks(1,1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);

	auto stats=new StatsWrapper(0.5,0.6,0.88,0.88);
	stats->add("c_{1}", simple->GetParameter("c_{1}"), 3, "mm #mus");
	stats->add("t_{0}", simple->GetParameter("t_{0}"), 4, "#mus");
//	stats->add("z_{offset}", offset, 3, "mm");
	//stats->addChiSquare(*simple);
	stats->draw();
	return means;
}
TH1* plotFittedTimeWalk( std::string filename="combinedFit.root", std::string object="quad/zResidualByToT") {
	auto uncorrected=getObjectFromFile<TH2>(object, filename);
	return plotFittedTimeWalk(uncorrected);
}
//draw four timewalk on 4 pads in one canvas
void combineFittedTimeWalkForChips(std::string filename="combinedFit.root", std::string object="zResidualByToT") {
	auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	canv->Divide(2,2);
	HistogramCombiner comb{"combinedTW"};
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto h=plotFittedTimeWalk(filename, chipDirectories[i]+"/"+object);
		//comb.add(h, "chip "+to_string(i+1) );
	}
	//comb.createCombined();
}


void combineTimeWalkResiduals(std::string filename="combinedFit.root", std::string object="zResidualByToTCorrected") {
	//auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	//canv->Divide(2,2);
	HistogramCombiner comb{"combinedTWResiduals"};
	for(int i=0;i<4;i++) {
		//canv->cd(i+1);
		auto corrected=getObjectFromFile<TH2>(chipDirectories[i]+"/"+object, filename);
		corrected->FitSlicesY();
		auto means=dynamic_cast<TH1*>( gDirectory->Get("zResidualByToTCorrected_1") );
		comb.add(means, "chip "+to_string(i+1) );
	}
	comb.createCombined();
}


//fitSlicesY with specification of range
TH1* fitDiffusionSlices(TH2* h2, std::string x="z") {
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	TF1* gausBG=new TF1("gaus","gaus(0)+[3]", -2,2);
	TF1* gausRange=new TF1("gausRange","gaus(0)", -0.5,0.5);
	gausRange->SetParameters(4E4,0.05,0.22);
//	TF1* exGaus=new TF1("exGaus", "[c]*[l]/2*exp([l]/2*(2*[m]+[l]*[s]*[s]-2*x))*TMath::Erfc( ([m]+[l]*[s]*[s]-x)/sqrt(2)/[s] )", -2, 2); //Exponentially modified gaussian distribution (see wiki)
//	exGaus->SetParameters(2E4,3.1,-0.25,0.2); // Constant, Lambda, Mean, Sigma
	//exGaus->FixParameter(1,3.1);
	gausBG->SetParameters(1E4,0,0.3,0); // Constant, Mean, Sigma, bg
	gausBG->SetParLimits(2, 1E-3,10);
	gaus->SetParameters(1E4,0,0.38); // Constant, Mean, Sigma, bg

	h2->FitSlicesY(x=="z" ? gausRange : gaus, 0/*firstbin*/, -1/*lastbin*/, 15/*min number of entries*/, "QNRM");
//	h2->FitSlicesY(gaus, 0, -1, 50, "QNR");

	//view background contributions
//	gDirectory->Get( (h2->GetName()+std::string("_3")).c_str() )->DrawClone();
//	new TCanvas();

	return dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
}


TF1* fitDiffusion( TH2* h2 , std::string x="x", double z0=-1, std::string canvname="canv") {
	double zmax=30;
	std::string LorT=x=="x" ? "T" : x=="z" ? "L" : x; //longitudinal or transverse

	TF1* drift=new TF1("drift", ("sqrt( pow([#sigma_{"+x+"0}],2) + pow([D_{"+x+"}],2)/10*(x-[z0]) )").c_str(), z0, zmax);
	drift->SetLineWidth(1);
	new TCanvas((canvname+"_"+x).c_str(), (canvname+"_"+x).c_str(), 800,600);

	fitDiffusionSlices(h2,x);

	TH1* h2_2=nullptr;
	if(false && x=="z") {
		TH1 *h22 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_1")).c_str() ));
		TH1 *h23 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_3")).c_str() ));
		h2_2=h23;//makeOperated(h22, h23, [](double l, double s) { return fabs(l) < 1E-30 ? 0 : sqrt( 1./(l*l)+s*s ); });
	} else {
		h2_2 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
	}
	if(!h2_2) { auto msg="could not get results from fit\n"; std::cerr<<msg; throw msg; }

	h2_2->GetXaxis()->SetTitle("z-position [mm]");
	h2_2->GetYaxis()->SetTitle( ("#sigma_{"+x+"} from fit to track-residual [mm]").c_str() );
	increaseAxisSize(h2_2, 0.05);	
	h2_2->GetYaxis()->SetRangeUser(0, x=="x"?0.8:0.8);
	h2_2->GetXaxis()->SetRangeUser(z0,zmax);
	
	//guess parameters
	if(x=="z") drift->FixParameter(2,z0);
	else drift->SetParameter(2,z0); //z0
	drift->SetParameter(1,0.3); //D
	if(x=="x") 	drift->FixParameter(0, 0.0158771);
	else drift->SetParameter(0, 0.15);//sigma0


	//add error because of guard
//	if(x=="x")
//	{
//		for(int i=12; i<=15; i++) {
//			h2_2->SetBinError(i, 0.01);
//		}
//	}

	//add error to fit
//	h2_2=addErrorToHist(h2_2, 1E-3); //set all error bins equal
	h2_2->Fit(drift, "", "", z0, zmax);

	gStyle->SetOptTitle(0);	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	h2_2->Draw();
	gPad->SetTicks(1,1);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	//nicer stat pane
	gStyle->SetOptFit(false);
	auto stats=new StatsWrapper();
	stats->add("D_{"+LorT+"}",  drift->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
	if(x!="x") stats->add("#sigma_{"+x+"0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" );
	stats->add("z0" ,  drift->GetParameter(2), 2, "mm" );
	stats->addChiSquare(*drift);
	stats->draw();
	
	h2_2->GetXaxis()->SetRangeUser(drift->GetParameter(2),zmax);

	//plot residuals
	const bool plotResiduals=false;
	if(plotResiduals) {
		auto residualHistogram=getResidualHistogram(h2_2, drift);
		new TCanvas(("residuals"+x).c_str(), "Residuals");
		residualHistogram->Draw();
	}
	
	return drift;

}


void plotDiffusionFromHist(std::string filename="combinedFit.root", std::string histogramName="locExp/xResidualByz", std::string dir="x") {
	auto hists=chipDirectories;
	hists.push_back("quad");
	for(std::string chip : hists )
//	std::string chip="quad";
	{
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histogramName, filename);
		h2->Rebin2D(4,1);
		fitDiffusion(h2, dir, -8, chip);
	}
}


void plotDiffusionUsingCut(
		std::string filename="combinedFit.root",
		std::string histogramName="locExp/xResidualByz",
		std::string dir="x", std::string alignFile="align.dat") {
//	auto hists=chipDirectories;
//	hists.push_back("quad");
//	for(std::string chip : hists )
	std::string chip="quad";
	Alignment align(alignFile);
	{
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histogramName, filename);
		h2=removeBinsByPosition(h2, [&](double x, double y) {
			return fabs(y) < 2.5*align.hitErrors.hitError(x).X();
		});
		new TCanvas();
		h2->Draw("colz");
		fitDiffusion(h2, dir, -0.8, chip);
	}
}

void plotDiffusionCombined(std::string filename="combinedFit.root") {
	TH2D* histogram=nullptr;
	for(std::string x : {"x", "z"}) {
		for(std::string chip : chipDirectories ) {
			TH2D* h=getObjectFromFile<TH2D>( chip+"/locExp"+x+"ResidualByz", filename);
			if(histogram) histogram->Add(h);
			else histogram=h;
		}
		fitDiffusion(histogram, x, 0);
	}
}

void plotZResidualsByToT(std::string filename="combinedFit.root", std::string histname="zResidualByToTCorrected") {
	HistogramCombiner combination("ToTCombination");
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	for(std::string chip : chipDirectories ) {
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histname, filename);

		h2->FitSlicesY(gaus, 0/*firstbin*/, -1/*lastbin*/, 5/*min number of entries*/, "QNR");

		auto h1=dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_1")).c_str() ));
		combination.add(h1, chip);
	}
	combination.createCombined();
}


TH1* plotSpot(std::string filename="fitted.root") {
	auto tree=getObjectFromFile<TTree>("fitResults",filename);
	new TCanvas("spot", "spot", 800,800);
	auto histForCounting=getHistFromTree(tree, "hitAverage.x", "hitAverage.x>0", "count", "");
	int nEntries=histForCounting->GetEntries();
	auto hist=getHistFromTree(tree, "residual.y/0.055:residual.x/0.055", "hitAverage.x>0 && flag>0", "spot(80,-40,40,80,-40,40)", "colz");
	hist->GetXaxis()->SetTitle("Columns");
	hist->GetYaxis()->SetTitle("Rows");
	hist->Scale(1./nEntries);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	return hist;
}


void plotSlicedDiffusionWithFit( std::string filename="combinedFit.root", double z0=-0.77, std::string object="quad/locExp/zResidualByzByToT" ) {
	auto hist3=getObjectFromFile<TH3D>(object, filename);
	auto zaxis=hist3->GetZaxis();


	//HistogramCombiner slices("slices");
	std::vector<std::string> slicename = {"0.15 #mus < ToT < 0.50 #mus", "ToT > 0.50 #mus"};
	std::vector<std::pair<double,double> > binRanges = { {7,21}, {22, zaxis->GetNbins()+1} }; //0.025 per bin!
	std::vector<TH1*> histograms= { nullptr, nullptr };
	auto stats=new StatsWrapper();
	TF1* drift=new TF1("drift", "sqrt( pow([#sigma_{z0}],2) + pow([D],2)/10*(x-[z0]) )", -1, 10);
	drift->SetLineWidth(1);
	TLegend* legend = new TLegend();
	for(int i : {0,1} ) {
		zaxis->SetRange(binRanges[i].first, binRanges[i].second);
		auto proj=(TH2*)hist3->Project3D("yx");

		auto hsigma=fitDiffusionSlices(proj, "z");

		auto name="slice "+std::to_string(i);
		histograms.at(i)=(TH1*)hsigma->Clone(name.c_str());

		histograms.at(i)->GetXaxis()->SetTitle("z-position [mm]");
		histograms.at(i)->GetYaxis()->SetTitle("#sigma_{z} from fit to track-residual [mm]" );
		//increaseAxisSize(histograms.at(i), 0.05);
		//histograms.at(i)->GetYaxis()->SetRangeUser(0,0.45);
		//histograms.at(i)->GetXaxis()->SetRangeUser(z0,23);

		legend->AddEntry(histograms.at(i), (slicename.at(i)+" ("+to_string_precision(proj->GetEntries()/hist3->GetEntries()*100,0)+"%)" ).c_str() );

		cout<<proj->GetEntries()<<" entries\n";

		//slices.add(clone,slicename.at(i)+";z-position [mm];#sigma_{z} from fit to track-residual [mm]");//+" ("+std::to_string(int(proj->GetEntries()/1000000+0.49))+"M hits)");

		if(i==1) {
			auto driftParams=fitDiffusion(proj, "z", z0);
			//cin.get();
			stats->add( "D_{L}",  driftParams->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
			stats->add( "#sigma_{z0}" ,  driftParams->GetParameter(0)*1E3, 0, "#mum" );
			stats->add( "z0" ,  driftParams->GetParameter(2), 2, "mm" );

			for(int j=0; j<4; j++) {
				drift->SetParameter(j, driftParams->GetParameter(j));
			}
		}
	}

	changeLegendStyle(legend, 1, 0.055);
	legend->SetX1NDC(0.18);

	//new TCanvas();
	//drift->Draw();
	histograms[0]->SetLineColor(kGreen+2);
	histograms[0]->Draw("same");
	histograms[1]->Draw("same");
	drift->Draw("same");
	legend->Draw();

	//slices.setStyle(10);slices.titleSize=0.05;
//	slices.setYRange({0,0.6});
//	slices.createCombined();

//	stats->draw();
	//gPad->SetMargin(0.15,0.1,0.15,0.1);

}

void combineFitDiffusion(
		std::string filename="combinedFit.root",
		std::vector<std::string> slicename = {"diffusion"},
		std::vector<std::string> objects = {"quad/locExp/xResidualByz"},
		double z0=-0.73 ) {

	std::vector<TH1*> histograms{objects.size(), nullptr};
	auto stats=new StatsWrapper();
	TF1* drift=new TF1("drift", "sqrt( pow([#sigma_{z0}],2) + pow([D],2)/10*(x-[z0]) )", -1, 10);
	drift->SetLineWidth(1);
	TLegend* legend = new TLegend();

	for(	int i=0; i<objects.size(); i++ ) {
		auto proj=getObjectFromFile<TH2D>(objects.at(i), filename);


		auto hsigma=fitDiffusionSlices(proj, "x");

		auto name="slice "+std::to_string(i);
		histograms.at(i)=(TH1*)hsigma->Clone(name.c_str());

		histograms.at(i)->GetXaxis()->SetTitle("z-position [mm]");
		histograms.at(i)->GetYaxis()->SetTitle("#sigma_{z} from fit to track-residual [mm]" );
		//increaseAxisSize(histograms.at(i), 0.05);
		//histograms.at(i)->GetYaxis()->SetRangeUser(0,0.45);
		//histograms.at(i)->GetXaxis()->SetRangeUser(z0,23);

		legend->AddEntry(histograms.at(i), (slicename.at(i)).c_str() );



	}

	if(objects.size()==2) {
		auto subtracted=makeOperated(histograms[0], histograms[1], [](double a,double b){return 2*b-a;});
		subtracted->Draw();
		histograms.push_back(subtracted);
		legend->AddEntry(subtracted, "Subtracted");

		drift->SetParameter(2,z0);
		drift->FixParameter(0, 0.0158771);
		subtracted->Fit(drift, "", "", z0, 9);

		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(1);

		subtracted->Draw();
		gPad->SetTicks(1,1);
		gPad->SetMargin(0.15,0.1,0.15,0.1);

		//nicer stat pane
		gStyle->SetOptFit(false);
		//cin.get();
		stats->add( "D_{L}",  drift->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
		stats->add( "#sigma_{z0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" );
		stats->add( "z0" ,  drift->GetParameter(2), 2, "mm" );
	}

	changeLegendStyle(legend, 1, 0.055);
	legend->SetX1NDC(0.18);

	//new TCanvas();
	//drift->Draw();
//	histograms[0]->SetLineColor(kGreen+2);
	stats->draw();
	for(int i=0; i<histograms.size(); i++) {
		histograms[i]->Draw("same");
	}
	drift->Draw("same");
	legend->Draw();

}


void plotDriftVelocity( std::string filename="fitted.root", std::string alignFile="align.dat") {
	Alignment align(alignFile);
	auto file=openFile(filename);
	align.updateDriftSpeed(*file);
}


void plotDriftVelocityNicely( std::string filename="fitted.root", std::string alignFile="../align.dat") {
	auto tree=getObjectFromFile<TTree>("fitResults",filename);
	Alignment	alignment(alignFile);
	std::cout<<std::to_string(alignment.driftSpeed.value)<<"\n";
	auto vDrift=getHistFromTree(
			tree, "(50-hitAverage.z)/"+std::to_string(alignment.driftSpeed.value)+":(50-laser.z)",
			"fabs(hitAverage.x)+fabs(hitAverage.y)>0 && nHitsPassed>1 && chip==0",
			"vDrift(81,-0.25,40.25)", "prof");
	vDrift->SetLineWidth(0);
	vDrift->SetMarkerStyle(7);

	auto driftRelation=new TF1("driftRelation", "x/[vdrift]+[intercept]");
	driftRelation->SetParameters(1,1);
	vDrift->Fit(driftRelation);

	vDrift->GetXaxis()->SetTitle("Laser position [mm]");
	vDrift->GetYaxis()->SetTitle("Average time of arrival [ns]");
	increaseAxisSize(vDrift,0.05);

	gPad->SetMargin(0.15,0.1,0.15,0.1);
	driftRelation->SetLineWidth(1);
	gPad->SetTicks(1,1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);

	auto stats=new StatsWrapper();
	stats->add("v_{drift}", driftRelation->GetParameter("vdrift")*1E3, 1, "#mum/ns");
	stats->draw();
}


void testTransformBackFunction() {
	Alignment align{"../align.dat"};
	std::string fileName="fitted.root";
	TFile file(fileName.c_str(), "READ");
	TTreeReader reader("fitResults", &file);
	TTreeReaderArray<double> hitpositionx{reader, "hits.position.x"}, hitpositiony{reader, "hits.position.y"};
	
	while(std::cin.get()!='q' and reader.Next() ) {
		for(int i=0; i<hitpositionx.GetSize(); i++) {
				cout<<hitpositionx[i]<<" "<<hitpositiony[i]<<" --> ";
				cout<<align.transformBack( {hitpositionx[i], hitpositiony[i], 0} )<<"\n";
		}
	}

}


TProfile2D* getCombinedDeformations(std::string histName="quad/locExp/xResidualByPosition",
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	auto prof=getObjectFromFile<TProfile2D>(histName, fileNames.front());
	for(int i=1; i<fileNames.size(); i++) {
		auto p=getObjectFromFile<TProfile2D>(histName, fileNames[i]);
		prof->Add(p);
	}

	return prof;
}

std::vector<std::vector<double>> fitDeformationCorrections(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition"
) {
	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues={
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ } };

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto canv=new TCanvas();
	canv->Divide(2,2);

	std::cout<<"{";
	for(int i=0; i<4; i++) {
		canv->cd(i+1);
		gPad->SetTicks(1,1);

//		auto xResidualByPosition=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(chipDirectories[i]+"/"+histName,fileNames);

		auto xResidualByx=xResidualByPosition->ProfileX();

		for(int j=0; j<correction->GetNpar(); j++) correction->SetParameter( j, estimatedValues[i][j] );

		xResidualByx->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);
		xResidualByx->SetTitle(";x-position [mm];x-residual [mm]");

		std::cout<<"{";
		for(int j=0; j<correction->GetNpar(); j++) {
			if(j==correction->GetNpar()-1) std::cout<<correction->GetParameter(j)<<( i!=3 ? "},\n" : "}}\n");
			else std::cout<<correction->GetParameter(j)<<", ";
			estimatedValues[i][j]=correction->GetParameter(j);
		}
	}

	return estimatedValues;
}



void fitDeformationCorrectionsPerSliceZ(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues={
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ } };

	const std::string parNames[13] = { "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "offset" };


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto xByXZCanvas=new TCanvas("xByXZ", "xByXZ", 1000,1000);
	xByXZCanvas->Divide(2,2);
	auto fitParametersCanvas=new TCanvas("fitSlices", "fitSlices", 1000,1000);
	fitParametersCanvas->Divide(2,2);

	std::cout<<"{";
	for(int iChip=0; iChip<4; iChip++) {
		xByXZCanvas->cd(iChip+1);
		gPad->SetTicks(1,1);
		auto xResidualByXZ=getCombinedDeformations(chipDirectories[iChip]+"/locExp/xResidualByXZ",fileNames);
		xResidualByXZ->Rebin2D(4,4);
		xResidualByXZ->SetAxisRange(-2,10,"Y");
		xResidualByXZ->Draw("colz");

		std::vector<TH1*> fittedParams;
		for(int j=0; j<13; j++) {
			auto h=newProjectionHistogram(xResidualByXZ, "_c"+std::to_string(iChip)+"f"+std::to_string(j), "x");
			h->SetTitle(("Chip "+to_string(iChip)+";z-position [mm]; parameter '"+parNames[j]+"'" ).c_str() );
			fittedParams.push_back( h );
		}

		fitParametersCanvas->cd(iChip+1);
		int iSlice=0;
		forEachSlice(xResidualByXZ, "x", [&](TH1* proj) {
			iSlice++;
			if(proj->GetEntries()<1000) return;
			for(int j=0; j<correction->GetNpar(); j++) correction->SetParameter( j, estimatedValues[iChip][j] );
			proj->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);

			for(int j=0; j<correction->GetNpar()-1; j++) {
				fittedParams[j]->SetBinContent(iSlice, correction->GetParameter(j));
				fittedParams[j]->SetBinError( iSlice, correction->GetParError(j));
			}

//			proj->DrawCopy();
//			gPad->Update();
//			std::cin.get();
			return;
		});

		HistogramCombiner comb("comb", {fittedParams[0], fittedParams[3], fittedParams[6], fittedParams[9]});
		comb.setStyle(7);
		comb.makeLegend=false;
		comb.createCombinedOnCanvas(gPad);
	}
}

std::vector<std::vector<double>> fitDeformationCorrectionsPerSliceY(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues=fitDeformationCorrections(fileNames);

	const std::string parNames[13] = { "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "offset" };


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto xByXZCanvas=new TCanvas("xByXY", "xByXY", 1000,1000);
	xByXZCanvas->Divide(2,2);
	auto fitParametersCanvas=new TCanvas("fitSlices", "fitSlices", 1000,1000);
	fitParametersCanvas->Divide(2,2);

	for(int iChip=0; iChip<4; iChip++) {
//		std::cout<<"\n\n";
		xByXZCanvas->cd(iChip+1);
		gPad->SetTicks(1,1);
		auto xResidualByXY=getCombinedDeformations(chipDirectories[iChip]+"/locExp/xResidualByPosition",fileNames);
		xResidualByXY->Rebin2D(4,16);
		//xResidualByXZ->SetAxisRange(-2,10,"Y");
		xResidualByXY->Draw("colz0");

		std::vector<TH1*> fittedParams;
		for(int j=0; j<13; j++) {
			auto h=newProjectionHistogram(xResidualByXY, "_c"+std::to_string(iChip)+"f"+std::to_string(j), "x");
			h->SetTitle(("Chip "+to_string(iChip)+";y-position [mm]; parameter '"+parNames[j]+"'" ).c_str() );
			fittedParams.push_back( h );
		}

		const int selectedParameter=0;
		fitParametersCanvas->cd(iChip+1);
		int iSlice=0;
		forEachSlice(xResidualByXY, "x", [&](TH1* proj) {
			iSlice++;
			if(proj->GetEntries()<100) return;
			for(int j=0; j<correction->GetNpar(); j++) {
				if(j%3==selectedParameter)
					correction->SetParameter( j, estimatedValues[iChip][j] );
				else
					correction->FixParameter( j, estimatedValues[iChip][j] );
			}
			proj->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);

			for(int j=0; j<correction->GetNpar()-1; j++) {
				auto parameterValue=correction->GetParameter(j);
				fittedParams[j]->SetBinContent(iSlice, parameterValue );
				fittedParams[j]->SetBinError( iSlice, correction->GetParError(j));
			}

//			proj->DrawCopy();
//			gPad->Update();
//			std::cin.get();
			return;
		});

		estimatedValues[iChip].resize(33);
		estimatedValues[iChip][32]=estimatedValues[iChip][12];
		for(int iPar=0; iPar<4; iPar++) {
			auto fitResult = fittedParams[iPar*3+selectedParameter]->Fit("pol4", "QSM");
			auto aCoefficient=fitResult->GetParams()[0];
			estimatedValues[iChip][iPar*3+selectedParameter]=aCoefficient;
			for(int iFit=1; iFit<5; iFit++) {
				estimatedValues[iChip][12+iPar+4*(iFit-1)]=fitResult->GetParams()[iFit]/aCoefficient;
//				std::cout<<iFit<<": "<<fitResult->GetParams()[iFit]/aCoefficient<<"\n";
			}
		}

		HistogramCombiner comb("comb", {fittedParams[0+selectedParameter], fittedParams[3+selectedParameter], fittedParams[6+selectedParameter], fittedParams[9+selectedParameter]});
		comb.setStyle(7);
		comb.makeLegend=false;
		comb.createCombinedOnCanvas(gPad);
	}

	return estimatedValues;
}


std::vector<std::vector<double>> fitDeformationCorrections2D(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat"
) {
	auto correction=new TF2("correction",
			"(1+[b0]*y+[c0]*y*y+[d0]*y*y*y+[e0]*y*y*y*y)*breitwigner(0)+(1+[b1]*y+[c1]*y*y+[d1]*y*y*y+[e1]*y*y*y*y)*breitwigner(3) "
			"+"
			"(1+[b2]*y+[c2]*y*y+[d2]*y*y*y+[e2]*y*y*y*y)*breitwigner(6)+(1+[b3]*y+[c3]*y*y+[d3]*y*y*y+[e3]*y*y*y*y)*breitwigner(9)+[offset]");
	correction->SetLineWidth(1);
	correction->SetNpx(1000);
	correction->SetNpy(1000);

	std::vector<std::vector<double> > estimatedValues=fitDeformationCorrectionsPerSliceY(fileNames);
//	{
//			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/   1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b,c,d,e*/ 0 /*offset*/ },
//			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/   1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ },
//			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ },
//			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ } };

//	auto fitted1DValues = fitDeformationCorrections(fileNames, histName);
//
//	for(int i=0; i<4; i++) {
//		for(int j=0; j<fitted1DValues[i].size()-1; j++) {
//			estimatedValues[i][j]=fitted1DValues[i][j];
//		}
//		estimatedValues[i][32]=fitted1DValues[i][12];
//	}

	Alignment align(alignFile);


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto canv=new TCanvas("2dCorrectinoFit","2D fit of deformations",1000,1200);
//	canv->Divide(2,2);
	canv->SetMargin(0.1,0.2,0.1,0.05);

	std::cout<<"{";
	for(int i=0; i<4; i++) {
//		canv->cd(i+1);
		gPad->SetTicks(1,1);

		auto corners=align.chips[i].getChipCorners();
		auto xminmax=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& a, const TVector3& b){ return a.x() < b.x(); });
		auto yminmax=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& a, const TVector3& b){ return a.y() < b.y(); });

//		auto xResidualByPosition=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(chipDirectories[i]+"/"+histName,fileNames);
		xResidualByPosition->Rebin2D(4,4);
//		xResidualByPosition->SetAxisRange( xminmax.first->x(),  xminmax.second->x() + (i==1 ? -0.5 : 0), "X");
//		xResidualByPosition->SetAxisRange( yminmax.first->y(),  yminmax.second->y(), "Y");
//		std::cout<<"\n"<<"yrange="<<yminmax.first->y()<<" - "<<yminmax.second->y()<<"\n";
		xResidualByPosition->SetTitle(";x-position [mm]; y-position [mm]; mean x-residual [mm]");
		xResidualByPosition->SetTitleOffset(1.2, "y");
		xResidualByPosition->SetTitleOffset(1.7, "z");
		xResidualByPosition->Draw(i ? "colz0same" : "colz0");

		std::cout<<"/*from slices, not fitted:*/";
		for(int j=0; j<correction->GetNpar(); j++) correction->FixParameter( j, estimatedValues[i][j] );

		const double margin[]={-0.1,-0.1,-0.05,0};
		correction->SetRange( xminmax.first->x()+margin[i], yminmax.first->y()+margin[i], xminmax.second->x()-margin[i], yminmax.second->y()-margin[i] );

		xResidualByPosition->Fit(correction, "QSR", "");//, chipRange[i].min, chipRange[i].max);
		//xResidualByPosition->SetTitle(";x-position [mm];x-residual [mm]");

		std::cout<<"{";
		for(int j=0; j<correction->GetNpar(); j++) {
			if(j==correction->GetNpar()-1) std::cout<<correction->GetParameter(j)<<( i!=3 ? "},\n" : "}}\n");
			else std::cout<<correction->GetParameter(j)<<", ";
			estimatedValues[i][j]=correction->GetParameter(j);
		}
	}

	align.drawChipEdges(false);

	return estimatedValues;
}


void compareDeformations(std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	std::string dir("quad");
//	for(const auto& dir :chipDirectories)
	{
	//	auto xResidualByPosition=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(dir+"/locExp/xResidualByPosition",fileNames );
		auto xResidualByx=xResidualByPosition->ProfileX();

	//	auto xResidualByPosition_cor=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_corrected", fileName );
		auto xResidualByPosition_cor=getCombinedDeformations(dir+"/corrected/xResidualByPosition",fileNames);
		auto xResidualByx_cor=xResidualByPosition_cor->ProfileX();

		HistogramCombiner deform("deform"+dir);
		deform.add(xResidualByx, "Before correction;x-position [mm]; x-residual [mm]");
		deform.add(xResidualByx_cor, "After correction");
		deform.createCombined();
	}
}


//In QUAD Frame
std::vector< std::vector<TVector3> > getChipAreaFromEdge(double dist=2, std::string alignFile="align.dat") {
	std::vector< std::vector<TVector3> > areas;
	Alignment align(alignFile);
	for(int i=0; i<4; i++) {
				auto corners=align.chips[i].getChipCorners();
				TVector3 meanPos=std::accumulate(corners.begin(), corners.end(), TVector3());
				meanPos*=1.0/corners.size();
				areas.push_back({});
				for(const auto& c : corners) {
					TVector3 unit=(c-meanPos).Unit();
					areas.back().emplace_back( c-dist*unit );
				}
	}
	return areas;
}

void drawAreaPolyLine(const std::vector<TVector3>& corners, Color_t color=kBlack) {
	TPolyLine l;
	for(auto& corner : corners) {
		l.SetNextPoint(corner.x(), corner.y());
	}
	l.SetNextPoint(corners[0].x(), corners[0].y());
	l.SetLineColor(color);
	l.DrawClone();
}

//do frequency on an unweighted profile, i.e. all entreis should have same weight
TH1D* getFrequencyInAreaFromProfile(TProfile2D* original, std::vector<std::vector<TVector3>> areas, double min=-0.1, double max=0.1, int nBins=80, double entryWeight=1.0) {
	TH1D* frequencyHist=new TH1D(
		(original->GetName()+std::string("freq")).c_str(),
		(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
		nBins, min, max);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if( original->GetBinEntries(bin) <= 0 ) continue;
			double binContent= original->GetBinContent(i,j); //possible rounding error here
			double binError=original->GetBinError(i,j);
//			frequencyHistAll->Fill(binContent);
			TVector3 binPoint{original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j),0};
			for(auto a : areas ) if( isInArea(binPoint, a) ) {
				frequencyHist->Fill(binContent);
				break;
			}
//				frequencyPull->Fill(binContent/error);
		}
	}
	return frequencyHist;
}
TH1D* getFrequencyInArea(TH2D* original, std::vector<std::vector<TVector3>> areas, double min=-0.1, double max=0.1, int nBins=80, double entryWeight=1.0) {
	TH1D* frequencyHist=new TH1D(
		(original->GetName()+std::string("freq")).c_str(),
		(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
		nBins, min, max);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if( fabs(original->GetBinContent(bin)) < 1E-20 ) continue;
			double binContent= original->GetBinContent(i,j); //possible rounding error here
			double binError=original->GetBinError(i,j);
//			frequencyHistAll->Fill(binContent);
			TVector3 binPoint{original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j),0};
			for(auto a : areas ) if( isInArea(binPoint, a) ) {
				frequencyHist->Fill(binContent);
				break;
			}
//				frequencyPull->Fill(binContent/error);
		}
	}
	return frequencyHist;
}

void drawDeformationsWithAreas(TProfile2D* prof, std::string alignFile) {
	new TCanvas( (prof->GetName()+std::string("deformcanv")).c_str(), "deformations", 900,1000)	;

	prof->SetTitle("");
	prof->GetYaxis()->SetTitleOffset(1.3);
	prof->GetXaxis()->SetTitleOffset(1.1);
	prof->GetZaxis()->SetTitleOffset(1.6);
	prof->Draw("colz0");

	gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal(alignFile);

	auto areas=getChipAreaFromEdge(3, alignFile);
	for(const auto& a : areas) drawAreaPolyLine(a);
	//new TCanvas();

	gStyle->SetOptStat(1);
	HistogramCombiner combined(prof->GetName()+std::string("freq"));
	auto freq=getFrequencyHistogramFromProfile(prof, -0.2, 0.2);
	combined.add(freq, "All bins (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	auto freqArea=getFrequencyInAreaFromProfile(prof, areas, -0.2, 0.2);
	combined.add(freqArea, "Bins inside selected area (RMS="+std::to_string( int( freqArea->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.createCombined();
}
void drawDeformationsWithAreas(TH2D* h2, std::string alignFile) {
	new TCanvas( (h2->GetName()+std::string("deformcanv")).c_str(), "deformations", 900,1000)	;

	h2->SetTitle("");
	h2->GetYaxis()->SetTitleOffset(1.3);
	h2->GetXaxis()->SetTitleOffset(1.1);
	h2->GetZaxis()->SetTitleOffset(1.6);
	h2->Draw("colz1");

	gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal(alignFile);

	auto areas=getChipAreaFromEdge(3, alignFile);
	for(const auto& a : areas) drawAreaPolyLine(a);

	//new TCanvas();
	gStyle->SetOptStat(1);
	HistogramCombiner combined(h2->GetName()+std::string("freq"));
	auto freq=getFrequencyHistogramFromProfile(h2);
	combined.add(freq, "All bins (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	auto freqArea=getFrequencyInArea(h2, areas);
	combined.add(freqArea, "Bins inside selected area (RMS="+std::to_string( int( freqArea->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.createCombined();
}

void drawDeformationsWithAreas(
		std::string histName="quad/locExp/xResidualByPosition",
		std::string fileName="fitted.root",
		std::string alignFile="align.dat") {
	auto prof=getObjectFromFile<TProfile2D>(histName, fileName);
	prof->SetMinimum(-0.2);
	prof->SetMaximum(0.2);
	prof->Rebin2D(8,4);
	drawDeformationsWithAreas(prof,alignFile);
//	prof->Draw();
}

void combineDeformations(
		std::string histName="quad/locExp/xResidualByPosition",
		std::vector<std::string> fileNames={"fitted.root"},
		std::string alignFile="align.dat") {

	auto prof=getCombinedDeformations(histName, fileNames);

	//rebin
	prof->Rebin2D(4,4);
	removeBinsWithFewerEntries(prof, 800);

	drawDeformationsWithAreas(prof, alignFile);
}

void combineDeformationFrequencies(
		std::vector<std::string> histNames={"quad/locExp/xResidualByPosition",  "quad/corrected/xResidualByPosition" },
		std::vector<std::string> histTitles={"Before correction", "After correction"},
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	HistogramCombiner combined("combFreq");
	auto histTitle=histTitles.begin();
	for(auto& histName : histNames) {
		auto prof=getCombinedDeformations(histName, fileNames);
		//rebin
		prof->Rebin2D(4,4);
		removeBinsWithFewerEntries(prof, 1000);
		auto hist=getFrequencyHistogramFromProfile(prof);
		combined.add(hist, *histTitle++ + "( RMS = "+to_string( int( hist->GetRMS()*1000) )+" #mum)");

	}
	combined.setNcolumn(1);
	combined.createCombined();
}

/*
TProfile2D* correctDeformation(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat") {
	auto parameters=fitDeformationCorrections(fileNames, histName);
	auto xResidualByPosition=getCombinedDeformations("quad/"+histName,fileNames);

	Alignment align(alignFile);
	auto corrected=makeAppliedByPosition(xResidualByPosition, [&](double xres, double x, double y) {
		auto chipNumber=align.getChipNumber({x,y,0});
		if(chipNumber) xres-=deformationCorrection(chipNumber-1, x, parameters);
		return xres;
	}, "corrected");

	corrected->SetMinimum(-0.1);
	corrected->SetMaximum(0.1);

	//rebin
	corrected->Rebin2D(4,4);
	removeBinsWithFewerEntries(corrected, 800);

	drawDeformationsWithAreas(corrected);
	return corrected;
}


TProfile2D* correctDeformation2D(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat") {
	auto parameters=fitDeformationCorrections2D(fileNames); //fitDeformationCorrectionsPerSliceY(fileNames);
	auto xResidualByPosition=getCombinedDeformations("quad/"+histName,fileNames);

	xResidualByPosition->Rebin2D(4,4);

	Alignment align(alignFile);
	auto corrected=makeAppliedByPosition(xResidualByPosition, [&](double xres, double x, double y) {
		auto chipNumber=align.getChipNumber({x,y,0});
		if(chipNumber) xres-=deformationCorrection2D(chipNumber-1, x, y, parameters);
		return xres;
	}, "corrected");

	corrected->SetMinimum(-0.1);
	corrected->SetMaximum(0.1);

	//rebin
//	corrected->Rebin2D(4,4);
	removeBinsWithFewerEntries(corrected, 800);

	drawDeformationsWithAreas(corrected);
	return corrected;
}
*/



void plotDeformationsXZ( std::string fileName, std::string objectName="quad/corrected/xResidualByXZ") {
	auto prof2d=getObjectFromFile<TProfile2D>(objectName, fileName);
	prof2d->Rebin2D(6,2);
	removeBinsWithFewerEntries(prof2d, 500);
	gStyle->SetPalette(kRainBow);
	prof2d->SetMinimum(-0.1);
	prof2d->SetMaximum(0.1);
	prof2d->Draw("colz0");
	gPad->SetMargin(0.10,0.15,0.1,0.05);
	gPad->SetTicks(1,1);
	increaseAxisSize(prof2d);
	prof2d->GetZaxis()->SetTitleSize(0.05);
	prof2d->GetZaxis()->SetLabelSize(0.05);
	prof2d->GetYaxis()->SetTitleOffset(1.05);

	HistogramCombiner combined{"combinedxByXZ"};
	auto freq=getFrequencyHistogramFromProfile(prof2d);
	combined.add(freq, "Mean x-residual per bin (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.setStyle(3);
	combined.createCombined();
}

void compareChipEdges(std::vector<std::string> alignFiles) {
	TH2D axis{"axisobj", ";x-position [mm];y-position [mm]",1000,0,29.04,1000,-14.55,25.05 };
	new TCanvas("chipEdges", "chipEdges", 850,1000)	;
	axis.DrawClone();
	const std::vector<Color_t> colors={kBlack,kRed,kBlue,kGreen+1};
	for(int i=0; i<alignFiles.size(); i++) {
		Alignment alignment(alignFiles[i]);
		alignment.drawChipEdges(false, colors[i%colors.size()]);
	}
	Alignment nominal(alignFiles[0]);
	std::cout<<"\n";
	for(int i=1; i<alignFiles.size(); i++){
		std::cout<<alignFiles[i]<<"\n";
		Alignment alignment(alignFiles[i]);
		for(int j=0; j<4; j++) {
			auto corners=alignment.chips[j].getChipCorners();
			auto nominalCorners=nominal.chips[j].getChipCorners();
			std::cout<<"Chip "<<j<<"\n"
			"First pad: "<<1E3*(corners[0].X()-nominalCorners[0].X())<<", "<<1E3*(corners[0].Y()-nominalCorners[0].Y())<<"\n"
			"Second pad: "<<1E3*(corners[1].X()-nominalCorners[1].X())<<", "<<1E3*(corners[1].Y()-nominalCorners[1].Y())<<"\n";
		}	
		std::cout<<"\n";
	}
	gStyle->SetOptStat(0);
}


//				{"run863/align.dat", 150},
//				{"run860/align.dat", 200},
//				{"run859/align.dat", 280},
//				{"run861/align.dat", 350},
//				{"run862/align.dat", 400} }
//{
//				{"run882/align.dat", 350},
//				{"run881/align.dat", 330},
//				{"run880/align.dat", 300},
//				{"run875/align.dat", 280},
//				{"run876/align.dat", 250},
//				{"run877/align.dat", 230},
//				{"run878/align.dat", 200},
//				{"run879/align.dat", 180} } ) {
//{
//				{"run889/align.dat", 100},
//				{"run890/align.dat", 150},
//				{"run891/align.dat", 200},
//				{"run892/align.dat", 250},
//				{"run893/align.dat", 300},
//				{"run894/align.dat", 350} }
//{"run906/align.dat", 200},
//{"run907/align.dat", 250},
//{"run908/align.dat", 300},
//{"run909/align.dat", 350} } ) {
//{run 918 to 922 = premixed?
//{run 924 to 928 = jansbottle
//{run 930 to 934 = T3K
//{run 935 to 939 = 10% isobutane //938 broken
//{run 940 to 948 = 10% isobutane more range 150,200....450,500,600
//{run 949 to 954 = T2K new CF4 150..400
TGraph* makeDriftVelocityGraph(
		std::vector<std::pair< std::string, int> > alignFiles = {
				{"run918/align.dat", 150},
				{"run919/align.dat", 200},
				{"run920/align.dat", 250},
				{"run921/align.dat", 300},
				{"run922/align.dat", 350} } ) {
	auto graph=new TGraph();
	for(auto p : alignFiles) {
		Alignment align(p.first);
		graph->SetPoint(graph->GetN(), p.second, 1E3*align.driftSpeed.value );
	}
	auto axisObj=new TH2D("axis", ";Drift field [V/cm];Drift velocity [#mum/ns]", 100,0,1000,100,0,100);
	axisObj->Draw();
	graph->SetMarkerStyle(kFullSquare);
	graph->Draw("LP");

	//T2K 0 ppm H20
	double E0[]={50, 55.7096, 62.0711, 69.1591, 77.0565, 85.8557, 95.6597, 106.583, 118.754, 132.315, 147.424, 164.258, 183.015, 203.914, 227.199, 253.144, 282.05, 314.258, 350.144, 390.127, 434.676, 484.313, 539.617, 601.237, 669.893, 746.388, 831.619, 926.583, 1032.39, 1150.28, 1281.63, 1427.99, 1591.05, 1772.73, 1975.16, 2200.71, 2452, 2732, 3044, 3391.58, 3778.87, 4210.38, 4691.17, 5226.86, 5823.73, 6488.75, 7229.71, 8055.28, 8975.12, 10000};
	double v0[]={2.04002, 2.30769, 2.60762, 2.94039, 3.3058, 3.70531, 4.13644, 4.59139, 5.06494, 5.54505, 6.01434, 6.46397, 6.87083, 7.21634, 7.49079, 7.68099, 7.77745, 7.7829, 7.69528, 7.52265, 7.27845, 6.97324, 6.62151, 6.2424, 5.84655, 5.44708, 5.06098, 4.69316, 4.35382, 4.04373, 3.77089, 3.52857, 3.31774, 3.13511, 2.98414, 2.85478, 2.75293, 2.67841, 2.63084, 2.6058, 2.61428, 2.65564, 2.72656, 2.83663, 2.96708, 3.13061, 3.30522, 3.53504, 3.76503, 4.04092};
	//10% iC4H10, 0 ppm H20
//  double v0[]={1.13426, 1.29086, 1.46898, 1.66206, 1.87347, 2.09407, 2.33156, 2.58205, 2.82872, 3.06816, 3.31184, 3.53749, 3.74917, 3.94401, 4.11816, 4.2715, 4.39709, 4.50143, 4.57951, 4.62988, 4.6597, 4.65654, 4.63386, 4.58217, 4.51047, 4.42354, 4.31817, 4.21079, 4.10252, 3.98931, 3.89045, 3.79583, 3.72062, 3.65847, 3.59705, 3.55317, 3.52304, 3.48572, 3.46525, 3.4603, 3.46522, 3.46371, 3.4889, 3.52826, 3.60014, 3.71298, 3.86396, 4.05255, 4.26794, 4.53487};
	for(auto& i : v0) i*=10; //cm -> mm
	auto gf2Graph=new TGraph(50, E0, v0);
	gf2Graph->Draw("LP");


//  double E[]={50, 66.081, 87.3341, 115.423, 152.545, 201.606, 266.447, 352.142, 465.398, 615.08, 812.903, 1074.35, 1419.88, 1876.54, 2480.08, 3277.73, 4331.91, 5725.14, 7566.47, 10000};
	//drift velocity 6000 ppm H20:
//  double v[]={0.54939, 0.74586, 1.02719, 1.44666, 2.08051, 3.01983, 4.27608, 5.56567, 6.35007, 6.27079, 5.53335, 4.60376, 3.81359, 3.25426, 2.90413, 2.74681, 2.80015, 3.06236, 3.52143, 4.16269};
  //drift velocity 2000 ppm H20
//  double v[]={1.03551, 1.17939, 1.34458, 1.53515, 1.75495, 2.00793, 2.30147, 2.63475, 3.01314, 3.43331, 3.89135, 4.38178, 4.88644, 5.39418, 5.86708, 6.29969, 6.65412, 6.92256, 7.08704, 7.1334, 7.08137, 6.91197, 6.6763, 6.3607, 6.00515, 5.63007, 5.23613, 4.86558, 4.50842, 4.18424, 3.88932, 3.63344, 3.40702, 3.21127, 3.05272, 2.91446, 2.80918, 2.72642, 2.66808, 2.64493, 2.65017, 2.68884, 2.75914, 2.85589, 2.98712, 3.14448, 3.33179, 3.54655, 3.7875, 4.06347};
//	drift velocity 1000 ppm H2O
  double v[]={1.37499, 1.57466, 1.79315, 2.04984, 2.33979, 2.66556, 3.03264, 3.44287, 3.88414, 4.35511, 4.85085, 5.35086, 5.83835, 6.2916, 6.68703, 7.01652, 7.25569, 7.3931, 7.42597, 7.36332, 7.19864, 6.96308, 6.65831, 6.3058, 5.92458, 5.52782, 5.14411, 4.77022, 4.42094, 4.10984, 3.82356, 3.57551, 3.35815, 3.17399, 3.01324, 2.88606, 2.77993, 2.69901, 2.64809, 2.62159, 2.63531, 2.66455, 2.72483, 2.83644, 2.98124, 3.13443, 3.32197, 3.53536, 3.77726, 4.05841};
  //drift velocity 10% argon at 2000 ppm H20, 50 ppm O2, 200 ppm N2
//	double v[]={0.6801, 0.77449, 0.88236, 1.00842, 1.15476, 1.32191, 1.5209, 1.74083, 1.98782, 2.25913, 2.53657, 2.82624, 3.11767, 3.39736, 3.65268, 3.89981, 4.10769, 4.29175, 4.44482, 4.55851, 4.64116, 4.68638, 4.6998, 4.67785, 4.62808, 4.5455, 4.4547, 4.34372, 4.23091, 4.11412, 4.0051, 3.90894, 3.8186, 3.74442, 3.67799, 3.62502, 3.58732, 3.55316, 3.52875, 3.51816, 3.50531, 3.50389, 3.52209, 3.5589, 3.62825, 3.73734, 3.88371, 4.07257, 4.29768, 4.56251};
	for(auto& i : v) i*=10; //cm -> mm
	auto gfGraph=new TGraph(50, E0, v);
	gfGraph->Draw("LP");

//	//10 000 ppm?
//	double E2[]={50, 55.7096, 62.0711, 69.1591, 77.0565, 85.8557, 95.6597, 106.583, 118.754, 132.315, 147.424, 164.258, 183.015, 203.914, 227.199, 253.144, 282.05, 314.258, 350.144, 390.127, 434.676, 484.313, 539.617, 601.237, 669.893, 746.388, 831.619, 926.583, 1032.39, 1150.28, 1281.63, 1427.99, 1591.05, 1772.73, 1975.16, 2200.71, 2452, 2732, 3044, 3391.58, 3778.87, 4210.38, 4691.17, 5226.86, 5823.73, 6488.75, 7229.71, 8055.28, 8975.12, 10000};
//	double v2[]={0.35907, 0.40116, 0.45148, 0.50546, 0.56689, 0.64289, 0.72332, 0.81287, 0.92589, 1.05285, 1.19929, 1.38075, 1.59089, 1.84457, 2.14159, 2.50569, 2.91964, 3.39137, 3.90911, 4.44992, 4.96804, 5.42065, 5.77045, 5.97494, 6.02044, 5.92467, 5.70562, 5.41022, 5.06565, 4.71996, 4.38335, 4.07543, 3.80605, 3.56861, 3.36595, 3.19188, 3.05432, 2.96176, 2.87972, 2.83484, 2.82422, 2.84967, 2.90525, 2.99882, 3.11716, 3.27441, 3.44601, 3.6502, 3.89077, 4.1448};
//	for(auto& i : v2) i*=10; //cm -> mm
	//2000 ppm H2O
//	double v2[]={1.03551, 1.17939, 1.34458, 1.53515, 1.75495, 2.00793, 2.30147, 2.63475, 3.01314, 3.43331, 3.89135, 4.38178, 4.88644, 5.39418, 5.86708, 6.29969, 6.65412, 6.92256, 7.08704, 7.1334, 7.08137, 6.91197, 6.6763, 6.3607, 6.00515, 5.63007, 5.23613, 4.86558, 4.50842, 4.18424, 3.88932, 3.63344, 3.40702, 3.21127, 3.05272, 2.91446, 2.80918, 2.72642, 2.66808, 2.64493, 2.65017, 2.68884, 2.75914, 2.85589, 2.98712, 3.14448, 3.33179, 3.54655, 3.7875, 4.06347};
//	for(auto& i : v2) i*=10; //cm -> mm
//	auto gf3Graph=new TGraph(50, E0, v2);
//	gf3Graph->Draw("LP");

	gStyle->SetOptStat(0);
	gPad->SetTicks(1,1);


	auto data=(new TText())->DrawText(100, 25, "Data");
	auto gas=(new TText())->DrawText(100,15, "T2K");
	auto gash2o=(new TText())->DrawText(100,5, "T2K + 1000 ppm H2O");
	for(auto& g : {data, gas, gash2o} ) g->SetTextFont(42);

	return graph;
}

void combineDriftVelocityGraphsIsobutane() {
	auto g1=makeDriftVelocityGraph({
				{"run935/align.dat", 150},
				{"run936/align.dat", 200},
				{"run937/align.dat", 250},
				{"run939/align.dat", 350}
	});
	auto g2=makeDriftVelocityGraph({
				{"run940/align.dat", 150},
				{"run941/align.dat", 200},
				{"run942/align.dat", 250},
				{"run943/align.dat", 300},
				{"run944/align.dat", 350},
				{"run945/align.dat", 400},
				{"run946/align.dat", 450},
				{"run947/align.dat", 500},
				{"run948/align.dat", 600}
	});
	g1->SetLineColor(kBlue);
	g1->SetMarkerColor(kBlue);
	g1->Draw("LP");
}

TGraph* makeDiffusionCoefficientGraph(
		std::vector<std::pair< std::string, int> > fitFiles = {
				{"run924/fitted.root", 150},
				{"run925/fitted.root", 200},
				{"run926/fitted.root", 250},
				{"run927/fitted.root", 300},
				{"run928/fitted.root", 350} },
		std::string histname="quad/locExp/zResidualByz"
) {
	auto graph=new TGraph();
	for(auto& f : fitFiles) {
		auto hist=getObjectFromFile<TH2D>(histname, f.first);
		auto diffusion=fitDiffusion(hist, "y");
		graph->SetPoint(graph->GetN(), f.second, 1E3*diffusion->GetParameter(1) );
	}
	return graph;
}

void combineDiffusionGraphs() {
	//jans bottle
	auto mixed= makeDiffusionCoefficientGraph({
		{"run924/fitted.root", 150},
		{"run925/fitted.root", 200},
		{"run926/fitted.root", 250},
		{"run927/fitted.root", 300},
		{"run928/fitted.root", 350} });
	auto premixed= makeDiffusionCoefficientGraph({
		{"run918/fitted.root", 150},
		{"run919/fitted.root", 200},
		{"run920/fitted.root", 250},
		{"run921/fitted.root", 300},
		{"run922/fitted.root", 350} });
	AllCombiner<TGraph> c("c");
	c.add(mixed, "Mixed;Drift field [V/cm];Longitudinal diffusion coefficient D_{L} [#mum/#sqrt{cm}]");
	c.add(premixed, "Premixed");

	c.createCombined();

}

void makeNHitsGraph(
	std::vector<std::pair< std::string, int> > fitFiles = {
			{"run910/fitted.root", 550},
			{"run911/fitted.root", 600},
			{"run912/fitted.root", 650},
			{"run913/fitted.root", 700},
			{"run914/fitted.root", 750},
			{"run915/fitted.root", 800},
			{"run916/fitted.root", 850},
			{"run917/fitted.root", 900}}
) {
	AllCombiner<TGraph> comb("comb");
	for(int i=1; i<=4; i++) {
		auto graph=new TGraph();
		for(auto p : fitFiles) {
			auto hist=getObjectFromFile<TH1>("chip"+std::to_string(i)+"/nHits", p.first);
			//auto mean=getMeanFromSimpleGausFit(*hist);
			graph->SetPoint(graph->GetN(), p.second, hist->GetMean() );
		}
		graph->Draw(i==1 ? "ALP" : "LPsame");
		comb.add(graph, "chip "+std::to_string(i)+";Threshold [e-];Number of hits");
	}
	comb.createCombined();
}

