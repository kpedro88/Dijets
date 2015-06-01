#define TreeClass_cxx

//custom headers
#include "TreeClass.h"
#include "KAsymBins.h"

//ROOT headers
#include <TFile.h>
#include <TH1.h>

//STL headers
#include <vector>
#include <cmath>

using namespace std;

void TreeClass::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L TreeClass.C
//      Root > TreeClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;
	
	//enable histo errors
	TH1::SetDefaultSumw2(kTRUE);
	
	//setup bins
	KAsymBins bins;
	if(!bins.parsed) return;
	//open output ROOT file
	//TFile* file = TFile::Open("tree_test_dijet/hist_QCD.root","RECREATE");
	TFile* file = TFile::Open("tree_dijet/hist_QCD.root","RECREATE");
	file->cd();
	//setup histos
	bins.MakeHistos();
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry % 10000 == 0) cout << "Skimming " << jentry << "/" << nentries << endl;
		
		//loop over jet type
		for(int ijt = 0; ijt < bins.jtype.size(); ++ijt){
			//jet variables
			double dijetpt, jet0eta, jet1eta, asymmetry;
			
			if(bins.jtype[ijt]==0){ //reco jet
				if(abs(RecoDijetDeltaPhi)<=2.7) continue;
				dijetpt = RecoDijetPt;
				jet0eta = RecoJet->at(0).Eta();
				jet1eta = RecoJet->at(1).Eta();
				asymmetry = RecoAsymmetry;
			}
			else if(bins.jtype[ijt]==1){ //gen jet
				if(abs(GenDijetDeltaPhi)<=2.7) continue;
				dijetpt = GenDijetPt;
				jet0eta = GenJet->at(0).Eta();
				jet1eta = GenJet->at(1).Eta();
				asymmetry = GenAsymmetry;
			}
		
			//loop over jet qty bins
			int ipt, iet;
			bool bpt, bet;
			bpt = bet = false;
			
			for(ipt = 0; ipt < bins.pt.size()-1; ++ipt){
				if(dijetpt>=bins.pt[ipt] && dijetpt<bins.pt[ipt+1]){
					bpt = true;
					break;
				}
			}
			if(!bpt) continue;

			for(iet = 0; iet < bins.eta.size()-1; ++iet){
				if(abs(jet0eta)>=bins.eta[iet] && abs(jet0eta)<bins.eta[iet+1] && abs(jet1eta)>=bins.eta[iet] && abs(jet1eta)<bins.eta[iet+1]){
					bet = true;
					break;
				}
			}
			if(!bet) continue;
			
			//loop over alpha types
			for(int iat = 0; iat < bins.atype.size(); ++iat){
				double alpha;
				if(bins.jtype[ijt]==0){ //reco jet
					if(bins.atype[iat]==0) alpha = RecoAlpha; //std
					else if(bins.atype[iat]==1) alpha = RecoApar; //parallel
					else if(bins.atype[iat]==2) alpha = RecoAperp; //perpendicular
				}
				else if(bins.jtype[ijt]==1){ //gen jet
					if(bins.atype[iat]==0) alpha = GenAlpha; //std
					else if(bins.atype[iat]==1) alpha = GenApar; //parallel
					else if(bins.atype[iat]==2) alpha = GenAperp; //perpendicular
				}
				
				//loop over alpha bins & fill histos
				for(int ial = 0; ial < bins.alpha.size()-1; ++ial){
					//inclusive binning
					if(alpha>=0.0 && alpha<=bins.alpha[ial+1]){
						bins.asym_incl[ijt][iat][ipt][iet][ial]->Fill(asymmetry,weight);
						bins.alpha_incl[ijt][iat][ipt][iet][ial]->Fill(alpha,weight);
					}
					
					//exclusive binning
					if(alpha>=bins.alpha[ial] && alpha<bins.alpha[ial+1]){
						bins.asym_excl[ijt][iat][ipt][iet][ial]->Fill(asymmetry,weight);
						bins.alpha_excl[ijt][iat][ipt][iet][ial]->Fill(alpha,weight);
					}	
				}
			}
		}
		
		// if (Cut(ientry) < 0) continue;
	}
	
	//finalize
	file->Write();
	file->Close();
}
