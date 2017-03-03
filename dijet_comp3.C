//custom headers
#include "KCode/KAsymBins.h"
#include "KCode/KTrend.h"

//STL headers
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

//--------------------------------------------------
//functions to initialize common drawing options
OptionMap* initGlobalOpt(){
	OptionMap* globalOpt = new OptionMap();
	globalOpt->Set<string>("prelim_text","Simulation (preliminary)");
	globalOpt->Set<string>("lumi_text","13 TeV");
	globalOpt->Set<bool>("checkerr",false);
	globalOpt->Set<double>("canvasH",475);
	globalOpt->Set<double>("sizeLeg",22);
	globalOpt->Set<double>("sizeSymb",0.1);
	return globalOpt;
}
OptionMap* initLocalOpt(){
	OptionMap* localOpt = new OptionMap();
	localOpt->Set<bool>("ratio",false);
	localOpt->Set<bool>("logy",false);
	return localOpt;
}

//creates samples, groups, plots
void dijet_comp3(bool print=false, string psuff="png", string pdir="plots", TFile* file=NULL){
	//open histogram file
	//if(!file) file = TFile::Open("/uscms_data/d3/pedrok/SUSY2015/crab/CMSSW_7_3_1_patch2/src/Dijets/tree_dijet/hist_QCD.root");
	if(!file) file = TFile::Open("/uscms_data/d3/pedrok/SUSY2015/crab/CMSSW_7_3_1_patch2/src/Dijets/tree_dijet/hist_QCD_test.root");
	
	//setup bins
	//KAsymBins bins;
	KAsymBins bins("input/input_bins_test.txt");
	if(!bins.parsed) return;
	//get histos from file
	bins.MakeHistos(file);
	
	Color_t colors[] = {kBlack, kBlue, kMagenta+2, kRed, kCyan+2, kMagenta, kOrange+7, kYellow+3};
	
	//global options
	OptionMap* globalOpt = new OptionMap();
	globalOpt->Set<string>("psuff",psuff);
	globalOpt->Set<string>("pdir",pdir);
	
	//axes: jtype, atype, pt, eta, alpha
	for(int ijt = 0; ijt < bins.jtype.size(); ++ijt){
		//quantity info (hardcoded, jet algos in different order)
		KJetQuantity::alg jtype = KJetQuantity::PUPPI;
		
		//for(int iat = 0; iat < bins.atype.size(); ++iat){
		for(int iat = 0; iat < 1; ++iat){		
			//todo: make alpha bin a separate dimension in AsymBins
			for(int iab = 0; iab < 2; ++iab){
				//for(int iet = 0; iet < bins.eta.size()-1; ++iet){
				for(int iet = 0; iet < 1; ++iet){
					//keep track of children
					vector<KTrend*> children_pt;
					
					for(int ipt = 0; ipt < bins.pt.size()-1; ++ipt){
						//keep track of children
						vector<KTrend*> children_alpha;
						
						for(int ial = 0; ial < bins.alpha.size()-1; ++ial){
							//histos
							TH1F* h_asym;
							vector<TH1F*> h_asym_split;
							
							//get histos
							if(iab==(int)KAlphaQuantity::Incl){
								h_asym = bins.asym_incl[ijt][iat][ipt][iet][ial];
								h_asym_split = bins.asym_split_incl[ijt][iat][ipt][iet][ial];
							}
							else if(iab==(int)KAlphaQuantity::Excl){
								h_asym = bins.asym_excl[ijt][iat][ipt][iet][ial];
								h_asym_split = bins.asym_split_excl[ijt][iat][ipt][iet][ial];
							}
							
							//skip empty histos
							if(h_asym->GetEntries()==0) continue;
							
							//option map for trend
							OptionMap* opt = new OptionMap();
							opt->Set<string>("name","asym");
							//quantities
							vector<string> quantity_names;
							quantity_names.push_back("Jet"); quantity_names.push_back("Alpha"); quantity_names.push_back("Pt"); quantity_names.push_back("Eta");
							opt->Set<vector<string> >("quantities",quantity_names);
							opt->Set<KJetQuantity::alg>("jtype",jtype);
							opt->Set<double>("etamin",bins.eta[iet]);
							opt->Set<double>("etamax",bins.eta[iet+1]);
							opt->Set<double>("ptmin",bins.pt[ipt]);
							opt->Set<double>("ptmax",bins.pt[ipt+1]);
							opt->Set<KAlphaQuantity::alph>("atype",(KAlphaQuantity::alph)iat);
							opt->Set<KAlphaQuantity::alphbin>("abin",(KAlphaQuantity::alphbin)iab);
							if(iab==(int)KAlphaQuantity::Incl){
								opt->Set<double>("amin",0.0);
								opt->Set<double>("amax",bins.alpha[ial+1]);
								opt->Set<double>("amean",bins.alpha[ial+1]);
								opt->Set<double>("ameanE",0.0);
								h_asym = bins.asym_incl[ijt][iat][ipt][iet][ial];
								h_asym_split = bins.asym_split_incl[ijt][iat][ipt][iet][ial];
							}
							else if(iab==(int)KAlphaQuantity::Excl){
								TH1F* h_alpha_excl = bins.alpha_excl[ijt][iat][ipt][iet][ial];
								opt->Set<double>("amin",bins.alpha[ial]);
								opt->Set<double>("amax",bins.alpha[ial+1]);
								opt->Set<double>("amean",h_alpha_excl->GetMean());
								opt->Set<double>("ameanE",h_alpha_excl->GetMeanError());
								h_asym = bins.asym_excl[ijt][iat][ipt][iet][ial];
								h_asym_split = bins.asym_split_excl[ijt][iat][ipt][iet][ial];
							}
							
							//option map for dist
							OptionMap* opt_d = new OptionMap();
							opt_d->Set<KDist::type>("type",KDist::Histo);
							opt_d->Set<KDist::source>("source",KDist::Object);
							opt_d->Set<TH1F*>("dist",h_asym);
							opt_d->Set<vector<TH1F*> >("hsplit",h_asym_split);
							opt_d->Set<vector<float> >("weights",bins.weights);
							opt_d->Set<string>("printname","asym");
							//hists already formatted, so formatting options not needed
							opt_d->Set<string>("xtitle","Asymmetry");
							//set plot local and global options
							OptionMap* plotGlobalOpt = initGlobalOpt();
							plotGlobalOpt->Set<int>("npanel",3);
							plotGlobalOpt->Set<double>("sizeLeg",20);
							plotGlobalOpt->Set<bool>("balance_panels",false);
							opt_d->Set<OptionMap*>("plotGlobalOpt",plotGlobalOpt);
							OptionMap* plotLocalOpt = initLocalOpt();
							plotLocalOpt->Set<bool>("logy",true);
							opt_d->Set<OptionMap*>("plotLocalOpt",plotLocalOpt);
							
							//option maps for fits
							vector<KNamed*> fitLines;
							OptionMap* opt_f = new OptionMap();
							opt_f->Set<Color_t>("color",kBlue);
							opt_f->Set<int>("line",2);
							opt_f->Set<string>("legname","Gaussian");
							opt_f->Set<bool>("do_ftail",false);
							KNamed* ftmp = new pair<string,OptionMap*>("gaussiancore",opt_f);
							fitLines.push_back(ftmp);
							
							opt_f = new OptionMap();
							opt_f->Set<Color_t>("color",kRed);
							opt_f->Set<int>("line",1);
							opt_f->Set<int>("panel",2);
							opt_f->Set<string>("legname","Crystal Ball");
							ftmp = new pair<string,OptionMap*>("crystalball",opt_f);
							fitLines.push_back(ftmp);
							
							//add fits to dist
							opt_d->Set<vector<KNamed*> >("fitLines",fitLines);
							
							//add dists to trend map
							vector<KNamed*> distLines;
							KNamed* dtmp = new pair<string,OptionMap*>("asym",opt_d);
							distLines.push_back(dtmp);
							opt->Set<vector<KNamed*> >("distLines",distLines);
							
							//make trend, draw self+fits
							KTrend* trend_asym = new KTrend(opt,globalOpt);
							trend_asym->Initialize();
							trend_asym->Draw(true,true,false,false,print);
							
							//store child
							children_alpha.push_back(trend_asym);
						} //loop over alpha
						/*
						//option map for trend
						//quantities inherited from children
						OptionMap* opt = new OptionMap();
						opt->Set<string>("name","extrap");
						opt->Set<vector<KTrend*> >("children",children_alpha);
						
						//dists from fit parameters of children, named dist:fit:fitqty
						//e.g. here, asym:gaussiancore:mu
						vector<KNamed*> distLines;
						string distpre = "asym:gaussiancore:";
						string distqty[] = {"mu","sigma","fmc","fdata","sf"};
						for(int d = 0; d < 5; ++d){
							//option map for dist
							OptionMap* opt_d = new OptionMap();
							opt_d->Set<KDist::type>("type",KDist::Graph);
							opt_d->Set<KDist::source>("source",KDist::Children);
							//formatting options
							opt_d->Set<Color_t>("color",kBlack);
							opt_d->Set<int>("marker",20);
							opt_d->Set<KLegend::Horz>("leg_horz",KLegend::left);
							opt_d->Set<KLegend::Vert>("leg_vert",KLegend::top);
							//fixed axis limits for sf
							if(distqty[d]=="sf") {
								opt_d->Set<bool>("yfixed",true);
								opt_d->Set<double>("ymin",0.5);
								opt_d->Set<double>("ymax",1.5);
							}
							//set plot local and global options
							OptionMap* plotGlobalOpt = initGlobalOpt();
							plotGlobalOpt->Set<int>("npanel",2);
							plotGlobalOpt->Set<bool>("balance_panels",false);
							opt_d->Set<OptionMap*>("plotGlobalOpt",plotGlobalOpt);
							OptionMap* plotLocalOpt = initLocalOpt();
							opt_d->Set<OptionMap*>("plotLocalOpt",plotLocalOpt);
							
							//each dit gets a LinOrConst fit applied
							vector<KNamed*> fitLines;
							OptionMap* opt_f = new OptionMap();
							//require const for sf
							if(distqty[d]=="sf") opt_f->Set<bool>("const",true);
							KNamed* ftmp = new pair<string,OptionMap*>("linorconst",opt_f);
							fitLines.push_back(ftmp);
							opt_d->Set<vector<KNamed*> >("fitLines",fitLines);
							
							//store dist with source qty
							opt_d->Set<string>("source_qty",distpre+distqty[d]);
							KNamed* dtmp = new pair<string,OptionMap*>("mc",opt_d);
							distLines.push_back(dtmp);
						}
						
						//add dists to trend map
						opt->Set<vector<KNamed*> >("distLines",distLines);
						
						//make trend, draw self+fits
						KTrend* trend_alpha = new KTrend(opt,globalOpt);
						trend_alpha->Initialize();
						trend_alpha->Draw(true,true,false,false,print);
						
						//store child
						children_pt.push_back(trend_alpha);
						*/
					} //loop over pt
					/*
					//option map for trend
					//quantities inherited from children
					OptionMap* opt = new OptionMap();
					opt->Set<string>("name","trend");
					opt->Set<vector<KTrend*> >("children",children_pt);
					
					//dists from fit parameters of children, named dist:fit:fitqty
					//e.g. here, extrap:linorconst:mu0
					vector<KNamed*> distLines;
					string distpre = "mc:linorconst:";
					string distqty[] = {"mu0","sigma0","fmc0","fdata0","sf0"};
					for(int d = 0; d < 5; ++d){
						//option map for dist
						OptionMap* opt_d = new OptionMap();
						opt_d->Set<KDist::type>("type",KDist::Graph);
						opt_d->Set<KDist::source>("source",KDist::Children);
						//formatting options
						opt_d->Set<Color_t>("color",kBlack);
						opt_d->Set<int>("marker",20);
						opt_d->Set<KLegend::Horz>("leg_horz",KLegend::right);
						opt_d->Set<KLegend::Vert>("leg_vert",KLegend::top);
						//fixed axis limits for sf
						if(distqty[d]=="sf0") {
							opt_d->Set<bool>("yfixed",true);
							opt_d->Set<double>("ymin",0.5);
							opt_d->Set<double>("ymax",1.5);
						}
						//set plot local and global options
						OptionMap* plotGlobalOpt = initGlobalOpt();
						plotGlobalOpt->Set<int>("npanel",2);
						plotGlobalOpt->Set<bool>("balance_panels",false);
						opt_d->Set<OptionMap*>("plotGlobalOpt",plotGlobalOpt);
						OptionMap* plotLocalOpt = initLocalOpt();
						opt_d->Set<OptionMap*>("plotLocalOpt",plotLocalOpt);
						
						//no fits here
						
						//store dist
						opt_d->Set<string>("source_qty",distpre+distqty[d]);
						KNamed* dtmp = new pair<string,OptionMap*>("mc",opt_d);
						distLines.push_back(dtmp);
					}
					//add dists to trend map
					opt->Set<vector<KNamed*> >("distLines",distLines);
					
					//make trend, draw self+fits
					KTrend* trend_pt = new KTrend(opt,globalOpt);
					trend_pt->Initialize();
					trend_pt->Draw(true,true,false,false,print);
					*/
				} //loop over eta
			} //loop over alpha binning
		} //loop over alpha type
	} //loop over jet type
	
	if(print){
		system(("zip -qr "+pdir+".zip "+pdir+"/").c_str());
	}
	
}
