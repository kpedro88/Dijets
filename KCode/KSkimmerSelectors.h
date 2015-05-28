#ifndef KSKIMMERSELECTORS_H
#define KSKIMMERSELECTORS_H

//custom headers
#include "Analysis/KCode/KSelection.h"
#include "KSkimmer.h"
#include "Analysis/KCode/KMath.h"

//ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

//STL headers
#include <string>
#include <vector>
#include <cstdlib>

//needed to write vector<TLorentzVector> to tree
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

using namespace std;

//base class for Selectors is in KSelection.h

//----------------------------------------------------
//stores event info
class KEventInfoSelector : public KSelector<KSkimmer> {
	public:
		//constructor
		KEventInfoSelector() : KSelector<KSkimmer>() { }
		KEventInfoSelector(string name_, OptionMap* localOpt_) : KSelector<KSkimmer>(name_,localOpt_) { canfail = false; }
		
		//accessors
		virtual void SetTree(TTree* tree_) { 
			tree = tree_;
			//set tree branches here
			tree->Branch("run",&run,"run/I");
			tree->Branch("lumi",&lumi,"lumi/I");
			tree->Branch("evt",&evt,"evt/I");
			tree->Branch("weight",&weight,"weight/F");
			//default values for variables
			run = lumi = evt = 0;
			weight = 0.0;
		}
		
		//set variable values for tree
		virtual bool Cut() {
			run = looper->run;
			lumi = looper->lumi;
			evt = looper->evt;
			weight = looper->weight;
			
			return true;
		}
		
		//member variables
		int run, lumi, evt;
		float weight;
};

//----------------------------------------------------
//creates jet 4vecs
class KJetSelector : public KSelector<KSkimmer> {
	public:
		//constructor
		KJetSelector() : KSelector<KSkimmer>() { }
		KJetSelector(string name_, OptionMap* localOpt_) : KSelector<KSkimmer>(name_,localOpt_) { 
			canfail = false;
			if(name.compare(0,3,"Gen")==0) do_gen = true;
			else do_gen = false;
		}
		
		//accessors
		virtual void SetTree(TTree* tree_) { 
			tree = tree_;
			string bpre = "";
			if(do_gen) bpre = "Gen";
			else bpre = "Reco";
			//set tree branches here
			tree->Branch((bpre+"Jet").c_str(),"std::vector<TLorentzVector>",&Jets,32000,0);
			//default values for variables
			Jets = NULL;
		}
		
		//used for non-dummy selectors
		virtual bool Cut() {
			//reset variables
			delete Jets; Jets = new vector<TLorentzVector>();
			
			//loop over jets
			unsigned jsize = 0;
			if(do_gen) {
				jsize = looper->refpt->size();
				Jets->reserve(jsize);
				for(unsigned j = 0; j < jsize; ++j){
					TLorentzVector tmp;
					tmp.SetPtEtaPhiE(looper->refpt->at(j),looper->refeta->at(j),looper->refphi->at(j),looper->refe->at(j));
					Jets->push_back(tmp);
				}
			}
			else {
				jsize = looper->jtpt->size();
				Jets->reserve(jsize);
				for(unsigned j = 0; j < jsize; ++j){
					TLorentzVector tmp;
					tmp.SetPtEtaPhiE(looper->jtpt->at(j),looper->jteta->at(j),looper->jtphi->at(j),looper->jte->at(j));
					Jets->push_back(tmp);
				}
			}
			
			return true;
		}
		
		//member variables
		bool do_gen;
		vector<TLorentzVector>* Jets;
};

//----------------------------------------------------
//selects events based on dijets
class KDijetSelector : public KSelector<KSkimmer> {
	public:
		//constructor
		KDijetSelector() : KSelector<KSkimmer>() { }
		KDijetSelector(string name_, OptionMap* localOpt_) : KSelector<KSkimmer>(name_,localOpt_) { 
			if(name.compare(0,3,"Gen")==0) do_gen = true;
			else do_gen = false;
		}
		
		//accessors
		virtual void SetTree(TTree* tree_) { 
			tree = tree_;
			string bpre = "";
			if(do_gen) bpre = "Gen";
			else bpre = "Reco";
			//set tree branches here
			tree->Branch((bpre+"Alpha").c_str(),&alpha,"alpha/D");
			tree->Branch((bpre+"Apar").c_str(),&apar,"apar/D");
			tree->Branch((bpre+"Aperp").c_str(),&aperp,"aperp/D");
			tree->Branch((bpre+"DijetPt").c_str(),&dijetpt,"dijetpt/D");
			tree->Branch((bpre+"DijetDeltaPhi").c_str(),&dijetdphi,"dijetdphi/D");
			tree->Branch((bpre+"Asymmetry").c_str(),&asymmetry,"asymmetry/D");
			//default values for variables
			alpha = apar = aperp = dijetpt = dijetdphi = asymmetry = 0.0;
		}
		virtual void SetSelection(KSelection<KSkimmer>* sel_) {
			sel = sel_;
			//set dependencies here
			if(do_gen){
				JetSel = sel->Get<KJetSelector*>("GenJet");
				if(!JetSel) {
					cout << "Input error: dependency GenJet failed in GenDijet!" << endl;
					depfailed = true;
				}
			}
			else {
				JetSel = sel->Get<KJetSelector*>("RecoJet");
				if(!JetSel) {
					cout << "Input error: dependency RecoJet failed in RecoDijet!" << endl;
					depfailed = true;
				}
			}
		}
		
		//used for non-dummy selectors
		virtual bool Cut() {
			//jet selector must be present for this to work
			if(depfailed) return false;
			
			//need 3 jets
			if(JetSel->Jets->size() < 3) return false;
			//minimum pT cut
			if(JetSel->Jets->at(2).Pt() < 10) return false;
			
			//alpha = extra jet activity: standard, parallel, perpendicular versions
			double dijet_sum_pt = JetSel->Jets->at(0).Pt() + JetSel->Jets->at(1).Pt();
			TLorentzVector dijet_axis = JetSel->Jets->at(0) - JetSel->Jets->at(1);
			TLorentzVector third_jet = JetSel->Jets->at(2);
			alpha = 2*third_jet.Pt()/dijet_sum_pt;
			apar = (2*third_jet*dijet_axis)/(dijet_axis.Pt()*dijet_sum_pt);
			aperp = sqrt(pow(alpha,2)-pow(apar,2));
			
			//dijet quantities: avg pt, delta phi, delta pt / sum pt
			dijetpt = dijet_sum_pt/2.;
			dijetdphi = KMath::DeltaPhi(JetSel->Jets->at(0).Phi(),JetSel->Jets->at(1).Phi());
			asymmetry = (JetSel->Jets->at(0).Pt() - JetSel->Jets->at(1).Pt())/dijet_sum_pt;
			
			return true;
		}
		
		//member variables
		bool do_gen;
		KJetSelector* JetSel;
		double alpha, apar, aperp, dijetpt, dijetdphi, asymmetry;
};

//-------------------------------------------------------------
//addition to KParser to create selectors
namespace KParser {
	template <>
	KSelector<KSkimmer>* processSelector<KSkimmer>(KNamed* tmp){
		KSelector<KSkimmer>* srtmp = 0;
		string sname = tmp->first;
		OptionMap* omap = tmp->second;
		
		//check for all known selectors
		if(sname=="EventInfo") srtmp = new KEventInfoSelector(sname,omap);
		if(sname=="GenJet" || sname=="RecoJet") srtmp = new KJetSelector(sname,omap);
		if(sname=="GenDijet" || sname=="RecoDijet") srtmp = new KDijetSelector(sname,omap);
		else {} //skip unknown selectors
		
		if(!srtmp) cout << "Input error: unknown selector " << sname << ". This selector will be skipped." << endl;
		
		return srtmp;
	}
}

#endif