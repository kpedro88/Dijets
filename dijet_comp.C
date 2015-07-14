//custom headers
#include "KCode/KAsymBins.h"
#include "KCode/KAsym.h"
#include "KCode/KDraw.h"

//STL headers
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

//creates samples, groups, plots
void dijet_comp(bool print=false, string psuff="png", string pdir="plots", TFile* file=NULL){
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
	
	//axes: jtype, atype, pt, eta, alpha
	for(int ijt = 0; ijt < bins.jtype.size(); ++ijt){
		for(int iat = 0; iat < bins.atype.size(); ++iat){
		//for(int iat = 0; iat < 1; ++iat){
			for(int iet = 0; iet < bins.eta.size()-1; ++iet){
			//for(int iet = 0; iet < 1; ++iet){
				//KAsymTrend here
				KAsymTrend* trend_excl = new KAsymTrend();
				KAsymTrend* trend_incl = new KAsymTrend();
				
				for(int ipt = 0; ipt < bins.pt.size()-1; ++ipt){
					//KAsymExtrap here
					KAsymExtrap* extrap_excl = new KAsymExtrap();
					KAsymExtrap* extrap_incl = new KAsymExtrap();
					
					for(int ial = 0; ial < bins.alpha.size()-1; ++ial){
						//get histos from file - exclusive alpha binning
						TH1F *h_asym_excl, *h_alpha_excl;
						vector<TH1F*> h_asym_split_excl;
						h_asym_excl = bins.asym_excl[ijt][iat][ipt][iet][ial];
						h_alpha_excl = bins.alpha_excl[ijt][iat][ipt][iet][ial];
						h_asym_split_excl = bins.asym_split_excl[ijt][iat][ipt][iet][ial];
						if(h_asym_excl->GetEntries()>0){ //skip histos with 0 entries						
							KAsymFit* asym_excl = new KAsymFit(h_asym_excl,h_asym_split_excl,bins.weights,(alg)bins.jtype[ijt],(alph)bins.atype[iat],bins.alpha[ial],bins.alpha[ial+1],h_alpha_excl->GetMean(),h_alpha_excl->GetMeanError(),bins.pt[ipt],bins.pt[ipt+1],bins.eta[iet],bins.eta[iet+1]);
							KDraw::DrawAsym(asym_excl,print,psuff,pdir);
							extrap_excl->push_back(asym_excl);
						}
						
						//get histos from file - inclusive alpha binning
						TH1F *h_asym_incl;
						vector<TH1F*> h_asym_split_incl;
						h_asym_incl = bins.asym_incl[ijt][iat][ipt][iet][ial];
						h_asym_split_incl = bins.asym_split_incl[ijt][iat][ipt][iet][ial];
						if(h_asym_incl->GetEntries()>0){ //skip histos with 0 entries
							KAsymFit* asym_incl = new KAsymFit(h_asym_incl,h_asym_split_incl,bins.weights,(alg)bins.jtype[ijt],(alph)bins.atype[iat],0.0,bins.alpha[ial+1],bins.alpha[ial+1],0.0,bins.pt[ipt],bins.pt[ipt+1],bins.eta[iet],bins.eta[iet+1]);
							KDraw::DrawAsym(asym_incl,print,psuff,pdir);
							extrap_incl->push_back(asym_incl);
						}
					}
					
					if(extrap_excl->asymfits.size()>1){
						KDraw::DrawExtrap(extrap_excl,bins.alpha[0],bins.alpha.back(),"excl",false,print,psuff,pdir);
						trend_excl->push_back(extrap_excl);
					}
					if(extrap_incl->asymfits.size()>1){
						KDraw::DrawExtrap(extrap_incl,bins.alpha[0],bins.alpha.back(),"incl",false,print,psuff,pdir);
						trend_incl->push_back(extrap_incl);
					}
				}
				
				KDraw::DrawExtrap(trend_excl,bins.pt[0],bins.pt.back(),"excl",true,print,psuff,pdir);
				KDraw::DrawExtrap(trend_incl,bins.pt[0],bins.pt.back(),"incl",true,print,psuff,pdir);
			}
		}
	}
	
	if(print){
		system(("zip -qr "+pdir+".zip "+pdir+"/").c_str());
	}
	
}
