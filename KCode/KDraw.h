#ifndef KDRAW_H
#define KDRAW_H

//custom headers
#include "KAsym.h"
#include "Analysis/KCode/KMap.h"
#include "Analysis/KCode/KLegend.h"
#include "Analysis/KCode/KPlot.h"

using namespace std;

namespace KDraw{

//--------------------------------------------------
//functions to initialize common drawing options
OptionMap* initGlobalOpt(){
	OptionMap* globalOpt = new OptionMap();
	globalOpt->Set<string>("prelim_text","Simulation (preliminary)");
	globalOpt->Set<string>("lumi_text","13 TeV");
	globalOpt->Set<bool>("checkerr",false);
	globalOpt->Set<double>("canvasH",475);
	return globalOpt;
}
OptionMap* initLocalOpt(){
	OptionMap* localOpt = new OptionMap();
	localOpt->Set<bool>("ratio",false);
	//localOpt->Set<bool>("logy",false);
	return localOpt;
}

//--------------------------------------------------
//function to draw asym histo with fit
void DrawAsym(KAsymFit* asym, bool print=false, string psuff="png", string pdir="plots"){
	//create base histo for drawing axes
	TH1F* hbase = (TH1F*)asym->hist->Clone("hbase");
	hbase->Reset();
	hbase->GetYaxis()->SetRangeUser(0.0,asym->hist->GetMaximum()*1.1);
	
	//make canvas/print name
	string oname = "asym__" + asym->printname;
	
	//get preamble text
	vector<string> preamble;
	asym->legnames[None] = "QCD";
	preamble.insert(preamble.begin(),asym->legnames.begin(),asym->legnames.end());
	if(preamble.back().size()>0) preamble.pop_back(); //ignore extra text if not used
	vector<int> extra_text_panels(preamble.size(),0);
	//add fit text
	preamble.insert(preamble.end(),asym->fitnames.begin(),asym->fitnames.end());
	vector<int> fit_panels(asym->fitnames.size(),1);
	extra_text_panels.insert(extra_text_panels.end(),fit_panels.begin(),fit_panels.end());
	
	//make plot options
	OptionMap* globalOpt = initGlobalOpt();
	globalOpt->Set<vector<string> >("extra_text",preamble);
	globalOpt->Set<int>("npanel",2);
	globalOpt->Set<bool>("balance_panels",0);
	globalOpt->Set<vector<int> >("extra_text_panels",extra_text_panels);
	globalOpt->Set<double>("sizeLeg",22);
	globalOpt->Set<double>("sizeSymb",0.05);
	OptionMap* localOpt = initLocalOpt();
	
	//make plot
	KPlot* plot = new KPlot(oname,localOpt,globalOpt);
	plot->Initialize(hbase);
	KLegend* kleg = plot->GetLegend();
	TCanvas* can = plot->GetCanvas();
	TPad* pad1 = plot->GetPad1();
	pad1->cd();

	//formatting
	asym->hist->SetLineWidth(2);
	kleg->AddHist(asym->hist);

	//build legend
	kleg->Build(KLegend::right);

	//draw blank histo for axes
	plot->DrawHist();
	
	//draw sample
	asym->hist->Draw("hist same");
	
	//linearized gaussian to show effect of tail
	TF1* gsn = new TF1("gsn","gaus",asym->hist->GetXaxis()->GetXmin(),asym->hist->GetXaxis()->GetXmax());
	gsn->SetParameters(asym->norm,asym->mu,asym->sigma);
	gsn->SetLineWidth(2);
	gsn->SetLineColor(kRed);
	gsn->SetLineStyle(2);
	gsn->Draw("same");
	
	//plot the classic gaussian fit (extended to the full range)
	TF1* gsnc = new TF1("gsnc","gaus",asym->hist->GetXaxis()->GetXmin(),asym->hist->GetXaxis()->GetXmax());
	gsnc->SetParameters(asym->gnorm,asym->gmu,asym->gsigma);
	gsnc->SetLineWidth(2);
	gsnc->SetLineColor(kBlue);
	gsnc->SetLineStyle(2);
	gsnc->Draw("same");
	
	//plot the novos fit as well
	asym->nfit->SetLineWidth(2);
	asym->nfit->SetLineColor(kRed);
	asym->nfit->SetLineStyle(1);
	asym->nfit->Draw("same");
	
	plot->GetHisto()->Draw("sameaxis"); //draw again so axes on top
	plot->DrawText();
	
	if(print){
		can->Print((pdir+"/"+oname+"."+psuff).c_str(),psuff.c_str());
	}
}


//--------------------------------------------------
//function to draw extrapolations or trends
void DrawExtrap(KAsymExtrap* extrap, double xmin, double xmax, string alphabin="", bool trend=false, bool print=false, string psuff="png", string pdir="plots"){
	//skip empty or otherwise pointless extraps
	if(!trend && extrap->asymfits.size()<2) return;
	
	//create graphs
	extrap->MakeGraphs(xmin,xmax);
	
	//loop over parameters to be extrapolated
	for(int p = 0; p < nExtrapPars; ++p){
		//make canvas/print name
		string oname;
		if(trend) oname = "trend_";
		else oname = "extrap_";
		if(alphabin.size()>0) oname += alphabin + "_";
		switch(p){
			case 0: oname += "mu"; break;
			case 1: oname += "sigma"; break;
			case 2: oname += "tau"; break;
			case 3: oname += "gmu"; break;
			case 4: oname += "gsigma"; break;
			case 5: oname += "fmc"; break;
			case 6: oname += "fdata"; break;
			case 7: oname += "sf"; break;
			default: oname += "";
		}
		oname += "_vs_";
		switch((qty)(extrap->q_varied)){
			case Jet: oname += "jet"; break;
			case Alpha:
				if((trend && extrap->asymextraps[0]->asymfits[0]->atype==Std) || extrap->asymfits[0]->atype==Std) { oname += "alpha"; }
				else if((trend && extrap->asymextraps[0]->asymfits[0]->atype==Par) || extrap->asymfits[0]->atype==Par) { oname += "apar"; }
				else if((trend && extrap->asymextraps[0]->asymfits[0]->atype==Perp) || extrap->asymfits[0]->atype==Perp) { oname += "aperp"; }
				break;
			case Pt: oname += "pt"; break;
			case Eta:  oname += "eta"; break;
			default:   oname += "";
		}
		string overall_common = extrap->GetCommonPrint();
		if(overall_common.size()>0) oname += "__" + overall_common;

		//create base histo for drawing axes
		//TH1F* hbase = new TH1F("hbase","",100,extrap->graph[p]->GetXaxis()->GetXmin(),extrap->graph[p]->GetXaxis()->GetXmax());
		TH1F* hbase = new TH1F("hbase","",100,xmin,xmax);
		hbase->GetXaxis()->SetTitle(extrap->graph[p]->GetXaxis()->GetTitle());
		hbase->GetYaxis()->SetTitle(extrap->graph[p]->GetYaxis()->GetTitle());
		double ymin = trend ? 0 : extrap->ymin[p];
		double yrange = extrap->ymax[p] - ymin;
		ymin = extrap->ymin[p]-yrange*0.1;
		if(trend && ymin < 0) ymin = 0;
		hbase->GetYaxis()->SetRangeUser(ymin,extrap->ymax[p]+yrange*0.1);
		//temporary
		if(p==7) hbase->GetYaxis()->SetRangeUser(0.5,1.5);

		//get preamble text - each time, b/c vector will be modified
		vector<string> preamble = extrap->GetCommonDescList();
		if(alphabin.size()>0) preamble.push_back(alphabin+" #alpha");
		vector<int> extra_text_panels(preamble.size(),0);
		//add fit text
		if(!trend){
			preamble.insert(preamble.end(),extrap->fitnames[p].begin(),extrap->fitnames[p].end());
			vector<int> fit_panels(extrap->fitnames[p].size(),1);
			extra_text_panels.insert(extra_text_panels.end(),fit_panels.begin(),fit_panels.end());
		}
		
		//make plot options
		OptionMap* globalOpt = initGlobalOpt();
		globalOpt->Set<vector<string> >("extra_text",preamble);
		if(!trend){
			globalOpt->Set<int>("npanel",2);
			globalOpt->Set<bool>("balance_panels",0);
			globalOpt->Set<vector<int> >("extra_text_panels",extra_text_panels);
		}
		globalOpt->Set<double>("sizeLeg",22);
		globalOpt->Set<double>("sizeSymb",0.05);
		OptionMap* localOpt = initLocalOpt();
		localOpt->Set<bool>("logy",false);
		//if(trend) localOpt->Set<bool>("logx",true);
		
		//make plot
		KPlot* plot = new KPlot(oname,localOpt,globalOpt);
		plot->Initialize(hbase);
		KLegend* kleg = plot->GetLegend();
		TCanvas* can = plot->GetCanvas();
		TPad* pad1 = plot->GetPad1();
		pad1->cd();
		
		//draw blank histo for axes
		plot->DrawHist();
		
		//draw extrap
		extrap->graph[p]->Draw("pe same");
		
		//draw fit
		if(!trend){
			extrap->gfit[p]->SetLineWidth(2);
			extrap->gfit[p]->SetLineColor(kRed);
			extrap->gfit[p]->SetLineStyle(1);
			extrap->gfit[p]->Draw("same");
		}
		
		//build legend
		kleg->AddHist(hbase); //for tick sizes
		if(trend) kleg->Build(KLegend::right,KLegend::top);
		else kleg->Build(KLegend::left,KLegend::top);

		//finish drawing
		plot->GetHisto()->Draw("sameaxis"); //draw again so axes on top
		plot->DrawText();
		
		if(print){
			can->Print((pdir+"/"+oname+"."+psuff).c_str(),psuff.c_str());
		}
	}
	
}

}

#endif