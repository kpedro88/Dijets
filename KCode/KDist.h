#ifndef KDIST_H
#define KDIST_H

//STL headers
#include <vector>
#include <string>
#include <sstream>

//ROOT headers
#include <TH1.h>

//custom headers
#include "Analysis/KCode/KMap.h"
#include "Analysis/KCode/KLegend.h"

using namespace std;

//forward declaration of fit class
class KFit;
class KFitQuantity;
typedef KMap<KFitQuantity*> FitQtyMap;

//---------------------------------------------------
//base class for KDist with common functions and vars
class KDist : public KQuantity {
	public:
		//enums
		enum type { NoType=0, Histo=1, Graph=2 };
		enum source { NoSource=0, Object=1, Tree=2, Children=3 };
		
		//constructors
		KDist() : KQuantity(), distSource(KDist::NoSource), distType(KDist::NoType), initialized(false) {}
		KDist(string name_, OptionMap* localOpt_) : KQuantity(name_,localOpt_), 
			distSource(KDist::NoSource), distType(KDist::NoType), source_qty(""), initialized(false),
			xtitle(""), ytitle(""),
			color(kBlack), marker(20), size(1.0), line(1), width(2), fit_quantities(new FitQtyMap())
		{
			//check source and type before initialization, in case the dist needs to be built
			localOpt->Get("source",distSource);
			localOpt->Get("type",distType);
			localOpt->Get("source_qty",source_qty);
		}
		//destructor
		virtual ~KDist() {}
		
		//common functions
		bool Initialized() { return initialized; }
		source GetSource() { return distSource; }
		string GetSourceQty() { return source_qty; }
		type GetType() { return distType; }
		FitQtyMap* GetFitQtys() { return fit_quantities; }
		vector<KFit*>& GetFits() { return fits; }

		//virtual functions (some will be specialized)
		virtual void Initialize() {}
		virtual TH1F* GetBaseHist() { return NULL; }
		virtual void AddToLegend(KLegend* kleg, string legname="") {}
		virtual void Draw(TPad* pad) {}
		
		//implemented in KFit.h to avoid circular dependency
		virtual void MakeFits(vector<KNamed*>& fitLines);
		
	protected:
		//member variables
		source distSource;
		type distType;
		string source_qty;
		bool initialized;
		string xtitle, ytitle;
		Color_t color;
		int marker;
		double size;
		int line, width;
		vector<KFit*> fits;
		FitQtyMap* fit_quantities;
};

//-----------------------------------------------
//class to store information about a distribution
//class T must be a ROOT histogram/graph/profile
template <class T>
class KDistT : public KDist {
	public:
		//constructor
		KDistT() : KDist(), dist(0) {}
		KDistT(string name_, OptionMap* localOpt_) : KDist(name_,localOpt_), dist(0)	{}
		//destructor
		virtual ~KDistT() {}
		
		//initialization from option map
		virtual void Initialize(){
			localOpt->Get("dist",dist);
			if(!dist) { //no distribution found
				initialized = false;
				return;
			}
			
			localOpt->Get("legname",legname);
			localOpt->Get("printname",printname);
			localOpt->Get("axisname",axisname);
			localOpt->Get("legvsname",legvsname);
			localOpt->Get("printvsname",printvsname);

			//apply settings to dist
			if(localOpt->Get("xtitle",xtitle)) dist->GetXaxis()->SetTitle(xtitle.c_str());
			if(localOpt->Get("ytitle",ytitle)) dist->GetYaxis()->SetTitle(ytitle.c_str());
			if(localOpt->Get("color",color)) { dist->SetMarkerColor(color); dist->SetLineColor(color); }
			if(localOpt->Get("marker",marker)) dist->SetMarkerStyle(marker);
			if(localOpt->Get("size",size))     dist->SetMarkerSize(size);
			if(localOpt->Get("line",line))     dist->SetLineStyle(line);
			if(localOpt->Get("width",width))   dist->SetLineWidth(width);
			
			//construct and apply fits
			vector<KNamed*> fitLines;
			bool hasFits = localOpt->Get("fitLines",fitLines);
			if(hasFits && fitLines.size()>0) MakeFits(fitLines);	
			
			initialized = true;
		}
		
		//accessors
		T* GetDist() { return dist; }
		void SetDist(T* dist_) { dist = dist_; }
		
	protected:
		//member variables
		T* dist;
};

//----------------------------------------------
//specific instance for histos
class KDistHisto : public KDistT<TH1F> {
	public:
		//constructor
		KDistHisto() : KDistT<TH1F>() {}
		KDistHisto(string name_, OptionMap* localOpt_) : KDistT<TH1F>(name_,localOpt_) {}
		//destructor
		virtual ~KDistHisto() {}
		
		//specialized functions
		TH1F* GetBaseHist(){
			TH1F* hbase = (TH1F*)dist->Clone("hbase");
			hbase->Reset();
			return hbase;
		}
		void AddToLegend(KLegend* kleg, string legname) { 
			kleg->AddEntry(dist,legname,"l");
		}
		void Draw(TPad* pad) {
			pad->cd();
			dist->Draw("hist same");
		}
};

//----------------------------------------------
//specific instance for graphs
class KDistGraph : public KDistT<TGraphAsymmErrors> {
	public:
		//constructor
		KDistGraph() : KDistT<TGraphAsymmErrors>() {}
		KDistGraph(string name_, OptionMap* localOpt_) : KDistT<TGraphAsymmErrors>(name_,localOpt_) {}
		//destructor
		virtual ~KDistGraph() {}
		
		//specialized functions
		TH1F* GetBaseHist(){
			TH1F* hbase = new TH1F("hbase","",1,dist->GetXaxis()->GetXmin(),dist->GetXaxis()->GetXmax());
			hbase->GetXaxis()->SetTitle(dist->GetXaxis()->GetTitle());
			hbase->GetYaxis()->SetTitle(dist->GetYaxis()->GetTitle());
			return hbase;
		}
		void AddToLegend(KLegend* kleg, string legname) { 
			kleg->AddEntry(dist,legname,"p"); //eventually "pe"
		}
		void Draw(TPad* pad) {
			pad->cd();
			dist->Draw("pe same");
		}
};

//addition to KParser to create dists
namespace KParser {
	KDist* processDist(KNamed* tmp){
		KDist* dtmp = 0;
		string name = tmp->first;
		OptionMap* omap = tmp->second;
		
		//check type
		KDist::type distType = KDist::NoType;
		omap->Get("type",distType);
		
		//construct if known type
		if(distType==KDist::Histo) dtmp = new KDistHisto(name,omap);
		else if(distType==KDist::Graph) dtmp = new KDistGraph(name,omap);
		
		if(!dtmp) cout << "Input error: unknown type " << distType << " for dist " << name << ". This dist will be skipped." << endl;
		
		return dtmp;
	}
}

#endif