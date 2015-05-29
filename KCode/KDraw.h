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
	
	//plot the fit as well
	asym->gfit->SetLineWidth(2);
	asym->gfit->SetLineColor(kRed);
	asym->gfit->SetLineStyle(1);
	asym->gfit->Draw("same");
	
	//linearized gaussian to show effect of cball tail
	TF1* gsn = new TF1("gsn","gaus",asym->hist->GetXaxis()->GetXmin(),asym->hist->GetXaxis()->GetXmax());
	gsn->SetParameters(asym->gfit->GetParameter(0),asym->mu,asym->sigma);
	gsn->SetLineWidth(2);
	gsn->SetLineColor(kRed);
	gsn->SetLineStyle(2);
	gsn->Draw("same");
	
	//line at start of tail
	Double_t bndR = asym->mu + asym->sigma*asym->a;
	TLine* line = new TLine(bndR,kleg->GetRange().first,bndR,kleg->GetRange().second);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->SetLineColor(kBlue);
	line->Draw("same");
	
	plot->GetHisto()->Draw("sameaxis"); //draw again so axes on top
	plot->DrawText();
	
	if(print){
		can->Print((pdir+"/"+oname+"."+psuff).c_str(),psuff.c_str());
	}
}

/*
//--------------------------------------------------
//function to draw resolution graphs
void DrawResolution(vector<KGroup*> groups, bool fit=true, bool print=false, string psuff="png", string pdir="plots"){
	//initial loop over groups
	KGroup* supergroup = new KGroup();
	for(int g = 0; g < groups.size(); g++){
		//create graph
		groups[g]->MakeGraph(fit);
		
		//check consistenty of varied qty
		if(g>0 && groups[g]->q_varied != groups[0]->q_varied){
			cout << "Error: inconsistent varied quantity among groups. DrawResolution() will exit now." << endl;
			return;
		}
		
		//create supergroup
		for(int s = 0; s < groups[g]->samples.size(); s++){
			supergroup->push_back(groups[g]->samples[s]);
		}
	}

	//check for overall commonalities
	for(int q = 0; q < supergroup->common.size(); q++){
		if(supergroup->common[q]){
			//remove overall commonality from individual groups
			//for labeling purposes
			for(int g = 0; g < groups.size(); g++){
				groups[g]->common[q] = false;
			}
		}
	}
	
	//make canvas/print name
	string oname = "reso";
	if(fit) oname += "_fit";
	oname += "_vs_";
	switch((qty)(groups[0]->q_varied)){
		case Algo: oname += "algo"; break;
		case Year: oname += "year"; break;
		case Lumi: oname += "lumi"; break;
		case Eta:  oname += "eta"; break;
		default:   oname += "";
	}
	oname += "__pt30";
	string overall_common = supergroup->GetCommonPrint();
	if(overall_common.size()>0) oname += "_" + overall_common;
	for(int g = 0; g < groups.size(); g++){
		oname += "__" + groups[g]->GetCommonPrint();
	}
	//todo: add "string extra" param to fn, to get eta bin info etc.?

	//create base histo for drawing axes
	TH1F* hbase = new TH1F("hbase","",100,groups[0]->graph->GetXaxis()->GetXmin(),groups[0]->graph->GetXaxis()->GetXmax());
	hbase->GetXaxis()->SetTitle(groups[0]->graph->GetXaxis()->GetTitle());
	hbase->GetYaxis()->SetTitle(groups[0]->graph->GetYaxis()->GetTitle());
	hbase->GetYaxis()->SetRangeUser(0.0,0.8);
	
	//get preamble text
	vector<string> preamble = supergroup->GetCommonDescList();
	//preamble.insert(preamble.begin(),"p_{T}^{Gen} > 10 GeV");
	preamble.insert(preamble.begin(),"#hat{p}_{T} = 30 GeV");

	//make plot options
	OptionMap* globalOpt = initGlobalOpt();
	globalOpt->Set<vector<string> >("extra_text",preamble);
	OptionMap* localOpt = initLocalOpt();
	
	//make plot
	KPlot* plot = new KPlot(oname,localOpt,globalOpt);
	plot->Initialize(hbase);
	KLegend* kleg = plot->GetLegend();
	TCanvas* can = plot->GetCanvas();
	TPad* pad1 = plot->GetPad1();
	pad1->cd();
	
	//draw blank histo for axes
	plot->DrawHist();
	
	//draw groups
	for(int g = 0; g < groups.size(); g++){
		groups[g]->graph->Draw("pe same");
		//add group to legend based on common cuts
		//kleg->AddEntry(groups[g]->graph,groups[g]->GetCommonDesc(),"pe"); //"e" option has wrong color until ROOT 5.34.11
		kleg->AddEntry(groups[g]->graph,groups[g]->GetCommonDesc(),"p");
	}

	//build legend
	kleg->AddHist(hbase); //for tick sizes
	kleg->Build(KLegend::left,KLegend::top);

	//finish drawing
	plot->GetHisto()->Draw("sameaxis"); //draw again so axes on top
	plot->DrawText();
	
	if(print){
		can->Print((pdir+"/"+oname+"."+psuff).c_str(),psuff.c_str());
	}
	
}
*/
}

#endif