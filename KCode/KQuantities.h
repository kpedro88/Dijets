#ifndef KQUANTITIES_H
#define KQUANTITIES_H

//STL headers
#include <iostream>

//custom headers
#include "KQuantity.h"

using namespace std;

//---------------------------------------------------
//version for alpha quantity (type, range, incl/excl)
class KAlphaQuantity : public KQuantity {
	public:
		//enums
		enum alph { Std=0, Par=1, Perp=2, alphsize=3 };
		enum alphbin { Incl=0, Excl=1, alphbinsize=2 };
	
		//constructor
		KAlphaQuantity() : KQuantity(), atype(Std), abin(Incl), amin(0), amax(0), amean(0), ameanE(0) { }
		KAlphaQuantity(string name_, OptionMap* localOpt_) : KQuantity(name_,localOpt_), atype(Std), abin(Incl), amin(0), amax(0), amean(0), ameanE(0) {
			//initialize from local opt
			localOpt->Get("atype",atype);
			localOpt->Get("abin",abin);
			localOpt->Get("amin",amin);
			localOpt->Get("amax",amax);
			localOpt->Get("amean",amean);
			localOpt->Get("ameanE",ameanE);
			
			//name options
			string aleg_opt[alphsize] = {"#alpha","#alpha_{#parallel}","#alpha_{#perp}"};
			string aleg = aleg_opt[(int)atype];
			string aprint_opt[alphsize] = {"alpha","apar","aperp"};
			string aprint = aprint_opt[(int)atype];
			string acl_opt[alphbinsize] = {"incl","excl"};
			string acl = acl_opt[(int)abin];
			
			//create names
			stringstream ss;
			ss << amin << " #leq " << aleg << " < " << amax;
			legname = ss.str();
			
			ss.str(string());
			ss << aprint << amin << "to" << amax;
			printname = ss.str();
			
			axisname = aleg + " (" + acl + ")";
			legvsname = acl + " " + aleg;
			printvsname = aprint + "_" + acl;
		}
		
		//accessors
		virtual double GetVal() { return amean; }
		virtual double GetErrLow() { return ameanE; }
		virtual double GetErrHigh() { return ameanE; }
		
		//comparison
		virtual bool Compare(const KQuantity& rhs){
			bool result = false;
			const KAlphaQuantity* p_rhs = dynamic_cast<const KAlphaQuantity*>(&rhs);
			if(p_rhs){
				result = (this->atype == p_rhs->atype && 
						  this->abin == p_rhs->abin && 
						  this->amin == p_rhs->amin && 
						  this->amax == p_rhs->amax);
			}
			return result;
		}
		
		//member variables
		alph atype;
		alphbin abin;
		double amin, amax;
		double amean, ameanE;
};

//---------------------------------------------------
//version for jet type
class KJetQuantity : public KQuantity {
	public:
		//enums
		enum alg { GEN=0, PUPPI=1, PFCHS=2, algsize=3 };
	
		//constructor
		KJetQuantity() : KQuantity(), jtype(GEN) { }
		KJetQuantity(string name_, OptionMap* localOpt_) : KQuantity(name_,localOpt_), jtype(GEN) {
			//initialize from local opt
			localOpt->Get("jtype",jtype);
			
			//create names
			string leg_opt[algsize] = {"GenJet","PUPPIJet","PFCHSJet"}; 
			legname = leg_opt[(int)jtype];
			string print_opt[algsize] = {"gen","puppi","pfchs"};
			printname = print_opt[(int)jtype];

			axisname = "Jet type";
			printvsname = "jet";
		}
		
		//accessors
		virtual double GetVal() { return (double)jtype; }
		virtual double GetErrLow() { return 0.0; }
		virtual double GetErrHigh() { return 0.0; }
		
		//comparison
		virtual bool Compare(const KQuantity& rhs){
			bool result = false;
			const KJetQuantity* p_rhs = dynamic_cast<const KJetQuantity*>(&rhs);
			if(p_rhs){
				result = (this->jtype == p_rhs->jtype);
			}
			return result;
		}
		
		//member variables
		alg jtype;
};

//---------------------------------------------------
//version for pT
class KPtQuantity : public KQuantity {
	public:
		//constructor
		KPtQuantity() : KQuantity(), ptmin(0), ptmax(0), pt(0) { }
		KPtQuantity(string name_, OptionMap* localOpt_) : KQuantity(name_,localOpt_), ptmin(0), ptmax(0), pt(0) {
			//initialize from local opt
			localOpt->Get("ptmin",ptmin);
			localOpt->Get("ptmax",ptmax);
			pt = (ptmin+ptmax)/2;
			
			//create names
			stringstream ss;
			ss << ptmin << " #leq p_{T} < " << ptmax;
			legname = ss.str();
			ss.str(string());
			ss << "pt" << ptmin << "to" << ptmax;
			printname = ss.str();

			axisname = "p_{T} [GeV]";
			printvsname = "pt";
		}
		
		//accessors
		virtual double GetVal() { return pt; }
		virtual double GetErrLow() { return pt-ptmin; }
		virtual double GetErrHigh() { return pt-ptmin; }
		
		//comparison
		virtual bool Compare(const KQuantity& rhs){
			bool result = false;
			const KPtQuantity* p_rhs = dynamic_cast<const KPtQuantity*>(&rhs);
			if(p_rhs){
				result = (this->ptmin == p_rhs->ptmin && this->ptmax == p_rhs->ptmax);
			}
			return result;
		}
		
		//member variables
		double ptmin, ptmax, pt;
};

//---------------------------------------------------
//version for eta
class KEtaQuantity : public KQuantity {
	public:
		//constructor
		KEtaQuantity() : KQuantity(), etamin(0), etamax(0), eta(0) { }
		KEtaQuantity(string name_, OptionMap* localOpt_) : KQuantity(name_,localOpt_), etamin(0), etamax(0), eta(0) {
			//initialize from local opt
			localOpt->Get("etamin",etamin);
			localOpt->Get("etamax",etamax);
			eta = (etamin+etamax)/2;
			
			//create names
			stringstream ss;
			ss << etamin << " #leq |#eta| < " << etamax;
			legname = ss.str();
			ss.str(string());
			ss << "eta" << etamin << "to" << etamax;
			printname = ss.str();

			axisname = "#eta";
			printvsname = "eta";
		}
		
		//accessors
		virtual double GetVal() { return eta; }
		virtual double GetErrLow() { return eta-etamin; }
		virtual double GetErrHigh() { return eta-etamin; }
		
		//comparison
		virtual bool Compare(const KQuantity& rhs){
			bool result = false;
			const KEtaQuantity* p_rhs = dynamic_cast<const KEtaQuantity*>(&rhs);
			if(p_rhs){
				result = (this->etamin == p_rhs->etamin && this->etamax == p_rhs->etamax);
			}
			return result;
		}
		
		//member variables
		double etamin, etamax, eta;
};

//-------------------------------------------------------------
//addition to KParser to create quantities
namespace KParser {
	KQuantity* processQuantity(KNamed* tmp){
		KQuantity* qty = 0;
		string name = tmp->first;
		OptionMap* omap = tmp->second;
		
		//check for all known quantities
		if(name=="Alpha") qty = new KAlphaQuantity(name,omap);
		else if(name=="Jet") qty = new KJetQuantity(name,omap);
		else if(name=="Pt") qty = new KPtQuantity(name,omap);
		else if(name=="Eta") qty = new KEtaQuantity(name,omap);
		else {} //skip unknown selectors
		
		if(!qty) cout << "Input error: unknown quantity " << name << ". This quantity will be skipped." << endl;
		
		return qty;
	}
}

#endif