#ifndef KFIT_H
#define KFIT_H

//STL headers
#include <vector>
#include <string>
#include <sstream>

//ROOT headers
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

//custom headers
#include "KQuantity.h"
#include "KDist.h"

using namespace std;

//---------------------------------------------------
//helper class to store fit results
class KFitQuantity : public KQuantity {
	public:
		//constructors
		KFitQuantity() : KQuantity(), value(0.), elow(0.), ehigh(0.), hasErr(true) {}
		KFitQuantity(string name_, string lname_, string pname_, string aname_, float value_, float elow_, float ehigh_) :
			KQuantity(name_,0), value(value_), elow(elow_), ehigh(ehigh_), hasErr(elow!=0 || ehigh!=0)
		{
			legname = lname_;
			printname = pname_;
			axisname = aname_;
		}
		//simpler constructor for symmetric errors
		KFitQuantity(string name_, string lname_, string pname_, string aname_, float value_, float err_) :
			KFitQuantity(name_,lname_,pname_,aname_,value_,err_,err_) {}
		//simpler constructor for case with no errors
		KFitQuantity(string name_, string lname_, string pname_, string aname_, float value_) :
			KFitQuantity(name_,lname_,pname_,aname_,value_,0.) {}
		
		//customized accessors
		string LegNameFull() {
			stringstream ss;
			int prec = 3;
			localOpt->Get("precision",prec);
			ss.precision(prec);
			ss << fixed << legname << " = " << value;
			if(hasErr) {
				if(elow==ehigh) ss << " #pm " << elow;
				else ss << "{}^{+" << ehigh << "}_{-" << elow << "}";
				//eventually use #plus and #minus?
			}
			return ss.str();
		}
		virtual double GetVal() { return value; }
		virtual double GetErrLow() { return elow; }
		virtual double GetErrHigh() { return ehigh; }
		bool HasError() { return hasErr; }
		
	protected:
		//member variables
		float value, elow, ehigh;
		bool hasErr;
};

typedef map<string,KFitQuantity*>::iterator FQMit;
typedef KMap<KFitQuantity*> FitQtyMap;

//---------------------------------------------------
//base class for KFit with common functions and vars
class KFit {
	public:
		//constructor
		KFit() : name(""), localOpt(new OptionMap()) {}
		KFit(string name_, OptionMap* localOpt_) : 
			name(name_), localOpt(localOpt_ ? localOpt_ : new OptionMap()),
			fit(NULL), printname(""), legname(""), color(kRed), line(1), width(2), fit_quantities(new FitQtyMap())
		{
			//initialization from option map
			localOpt->Get("printname",printname);
			localOpt->Get("legname",legname);
			localOpt->Get("color",color);
			localOpt->Get("line",line);
			localOpt->Get("width",width);
		}
		//destructor
		virtual ~KFit() {}
		
		//accessors
		FitQtyMap* GetFitQtys() { return fit_quantities; }
		void AddToLegend(KLegend* kleg, bool addToPrev){
			vector<string> extra_text;
			extra_text.reserve(fit_quantities->GetTable().size());
			for(unsigned q = 0; q < fit_quantities_list.size(); ++q){
				extra_text.push_back(fit_quantities_list[q]->LegNameFull());
			}
			//default is panel 1
			int panel = 1;
			localOpt->Get("panel",panel);
			kleg->AddEntry(fit,legname,"l",panel,extra_text,addToPrev);
		}
		virtual void Draw(TPad* pad) {
			//format fit
			fit->SetLineColor(color);
			fit->SetLineStyle(line);
			fit->SetLineWidth(width);
			
			pad->cd();
			fit->Draw("same");
		}
		
		//virtual functions
		virtual void SetDist(KDist* dist_) { }
		virtual KDist* GetDist() { return NULL; }
		virtual void Fit() {}
		
	protected:
		//member variables
		string name;
		OptionMap* localOpt;
		TF1* fit;
		string printname, legname;
		Color_t color;
		int line, width;
		vector<KFitQuantity*> fit_quantities_list;
		FitQtyMap* fit_quantities;
		
		//keep track of fit qtys for legend (in order) and searching (in map)
		void StoreFitQty(KFitQuantity* fq){
			//prepend fit name to qtys in map
			fit_quantities->Add(name+":"+fq->GetName(),fq);
			fit_quantities_list.push_back(fq);
		}
};

//---------------------------------------------------
//fit class template to align with KDist templates
//specific fits will inherit from a version of this
template <class T>
class KFitT : public KFit {
	public:
		//constructors
		KFitT() : KFit(), dist(NULL) {}
		KFitT(string name_, OptionMap* localOpt_) : KFit(name_,localOpt_), dist(NULL) {}
		//destructor
		virtual ~KFitT() {}
		
		//accessors
		virtual void SetDist(KDist* dist_) { dist = static_cast<KDistT<T>*>(dist_); }
		virtual KDist* GetDist() { return static_cast<KDist*>(dist); }
		
		//special "smart fit" procedure from JECs
		virtual TFitResultPtr SmartFit(T* dist_, TF1* fit_, int nrep=30){
			vector<double> params, errors;
			params.resize(fit_->GetNpar(),0);
			errors.resize(fit_->GetNpar(),0);
			double chi2tmp, ndftmp;
			TFitResultPtr result;
			
			double chi2best = -1;
			for(int f = 0; f < nrep; ++f){
				TFitResultPtr restmp = dist_->Fit(fit_,"NQRS");
				double ctmp = fit_->GetChisquare()/fit_->GetNDF();
				
				if(ctmp < chi2best || chi2best==-1){
					result = restmp;
					chi2best = ctmp;
					//get values from fit
					for(int p = 0; p < fit_->GetNpar(); ++p){
						params[p] = fit_->GetParameter(p);
						errors[p] = fit_->GetParError(p);
					}
					chi2tmp = fit_->GetChisquare();
					ndftmp = fit_->GetNDF();
				}
			}
			
			//put values back in fit
			for(int p = 0; p < fit_->GetNpar(); ++p){
				fit_->SetParameter(p,params[p]);
				fit_->SetParError(p,errors[p]);
			}
			fit_->SetChisquare(chi2tmp);
			fit_->SetNDF(ndftmp);
			
			//return fit result
			return result;
		}
		
	protected:
		//member variables
		KDistT<T>* dist;
};

//typedefs
typedef KFitT<TH1F> KFitHisto;
typedef KFitT<TGraphAsymmErrors> KFitGraph;

//specializations

//forward declaration of fit creator in KParser
namespace KParser {
	KFit* processFit(KNamed* tmp);
}

//implementation of KDist functions which depend on KFit class
//(to avoid circular dependency)
void KDist::MakeFits(vector<KNamed*>& fitLines){
	fits.reserve(fitLines.size());
	for(unsigned f = 0; f < fitLines.size(); ++f){
		KFit* ftmp = KParser::processFit(fitLines[f]);
		if(!ftmp) continue;
		
		ftmp->SetDist(this);
		
		//do fit and store qtys in dist
		//prepend dist name to qtys
		ftmp->Fit();
		FitQtyMap* fqm = ftmp->GetFitQtys();
		for(FQMit it = fqm->GetTable().begin(); it != fqm->GetTable().end(); ++it){
			fit_quantities->Add(name+":"+it->first,it->second);
		}

		//store fit
		fits.push_back(ftmp);
	}
}

#endif