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
	if(!file) file = TFile::Open("/uscms_data/d3/pedrok/SUSY2015/crab/CMSSW_7_3_1_patch2/src/Dijets/tree_dijet/hist_QCD.root");
	
	//setup bins
	KAsymBins bins;
	if(!bins.parsed) return;
	//setup histo names
	bins.MakeHistos(true);
	
	Color_t colors[] = {kBlack, kBlue, kMagenta+2, kRed, kCyan+2, kMagenta, kOrange+7, kYellow+3};
	
	//axes: jtype, atype, pt, eta, alpha
	for(int ijt = 0; ijt < bins.jtype.size(); ++ijt){
		for(int iat = 0; iat < bins.atype.size(); ++iat){
			for(int ipt = 0; ipt < bins.pt.size()-1; ++ipt){
				//KAsymTrend here?
				for(int iet = 0; iet < bins.eta.size()-1; ++iet){
					//KAsymExtrap here
					KAsymExtrap* extrap_excl = new KAsymExtrap();
					KAsymExtrap* extrap_incl = new KAsymExtrap();
					
					for(int ial = 0; ial < bins.alpha.size()-1; ++ial){
						//get histos from file - exclusive alpha binning
						TH1F *h_asym_excl, *h_alpha_excl;					
						h_asym_excl = (TH1F*)file->Get(bins.asym_excl_names[ijt][iat][ipt][iet][ial].c_str());
						h_alpha_excl = (TH1F*)file->Get(bins.alpha_excl_names[ijt][iat][ipt][iet][ial].c_str());
						if(h_asym_excl->GetEntries()>0){ //skip histos with 0 entries						
							KAsymFit* asym_excl = new KAsymFit(h_asym_excl,h_alpha_excl,(alg)bins.jtype[ijt],(alph)bins.atype[iat],bins.alpha[ial],bins.alpha[ial+1],bins.pt[ipt],bins.pt[ipt+1],bins.eta[iet],bins.eta[iet+1]);
							KDraw::DrawAsym(asym_excl,print,psuff,pdir);
							extrap_excl->push_back(asym_excl);

							//inclusive = exclusive for last bin
							if(ial==bins.alpha.size()-2) {
								extrap_incl->push_back(asym_excl);
								continue;
							}
						}
						
						//get histos from file - inclusive alpha binning
						TH1F *h_asym_incl, *h_alpha_incl;
						h_asym_incl = (TH1F*)file->Get(bins.asym_incl_names[ijt][iat][ipt][iet][ial].c_str());
						h_alpha_incl = (TH1F*)file->Get(bins.alpha_incl_names[ijt][iat][ipt][iet][ial].c_str());
						if(h_asym_incl->GetEntries()>0){ //skip histos with 0 entries
							KAsymFit* asym_incl = new KAsymFit(h_asym_incl,h_alpha_incl,(alg)bins.jtype[ijt],(alph)bins.atype[iat],bins.alpha[ial],bins.alpha.back(),bins.pt[ipt],bins.pt[ipt+1],bins.eta[iet],bins.eta[iet+1]);
							KDraw::DrawAsym(asym_incl,print,psuff,pdir);
							extrap_incl->push_back(asym_incl);
						}
					}
					
					KDraw::DrawExtrap(extrap_excl,"excl",print,psuff,pdir);
					KDraw::DrawExtrap(extrap_incl,"incl",print,psuff,pdir);
				}
			}
		}
	}
	
	if(print){
		system(("zip -qr "+pdir+".zip "+pdir+"/").c_str());
	}
	
}
