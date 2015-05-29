#ifndef KSAMPLE_H
#define KSAMPLE_H

//ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include "Math/MinimizerOptions.h"
#include <TVirtualFitter.h>

//STL headers
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace TMath;

enum qty { None=0, Jet=1, Alpha=2, Pt=3, Eta=4, Extra=5, qtySize=6 };
enum alg { Reco=0, Gen=1 };
enum alph { Std=0, Par=1, Perp=2 };

//------------------------------------------
//Right-sided Crystal Ball function
//parameters:
//N, mu, sigma, a, n
//0,  1,     2, 3, 4
Double_t cball(Double_t *x, Double_t *par){
	//ensure sigma > 0 and a > 0
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	par[3] = Abs(par[3]);
	Double_t a = par[3];
	par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
	Double_t n = par[4];
	Double_t arg = (x[0]-mu)/sigma;
	
	double d = n/a;

	//right tail
	if(arg >= a){
		return N*Exp(-Power(a,2)/2)*Power((d-a+arg)/d,-n);
	}
	//core
	else{
		return N*Exp(-Power(arg,2)/2);
	}

}

//------------------------------------------
//Right-sided log(Crystal Ball) function
//parameters:
//N, mu, sigma, a, n
//0,  1,     2, 3, 4
Double_t logcball(Double_t *x, Double_t *par){
	//ensure sigma > 0 and a > 0
	par[0] = Abs(par[0]);
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	par[3] = Abs(par[3]);
	Double_t a = par[3];
	par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
	Double_t n = par[4];
	Double_t arg = (x[0]-mu)/sigma;
	
	double d = n/a;

	//right tail
	if(arg >= a){
		return Log(N) - Power(a,2)/2 - n*Log((d-a+arg)/d);
	}
	//core
	else{
		return Log(N) - Power(arg,2)/2;
	}

}

//-----------------------------------------------------------------------------------------
//class to keep track of set of cut values, corresponding histo, fit, various string labels
class KAsymFit {
	public:
		//constructor
		KAsymFit(TH1F* hist_, TH1F* h_alpha, alg jtype_, alph atype_, double amin_, double amax_, double ptmin_, double ptmax_, double etamin_, double etamax_, string extra_="", Color_t color_=kBlack) :
			jtype(jtype_), atype(atype_), amin(amin_), amax(amax_), ptmin(ptmin_), ptmax(ptmax_), pt((ptmin+ptmax)/2), etamin(etamin_), etamax(etamax_), eta((etamin+etamax)/2), extra(extra_), color(color_),
			legnames(qtySize,""), printnames(qtySize,""), cutname(""), printname(""),
			hist(hist_), gfit(NULL), mean(0), meanE(0), rms(0), rmsE(0), Nevents(0), mu(0), muE(0), sigma(0), sigmaE(0), a(0), aE(0), n(0), nE(0), chi2ndf(0)
		{
			//create legend names (descriptions) & print names
			if(jtype==Reco) { legnames[Jet] = "RecoJet"; printnames[Jet] = "Reco"; }
			else if(jtype==Gen) { legnames[Jet] = "GenJet"; printnames[Jet] = "Gen"; }
			printname += printnames[Jet];
			
			stringstream ss;
			string aleg, aprint;
			if(atype==Std) {       aleg = "#alpha";             aprint = "Alpha"; }
			else if(atype==Par) {  aleg = "#alpha_{#parallel}"; aprint = "Apar"; }
			else if(atype==Perp) { aleg = "#alpha_{#perp}";     aprint = "Aperp"; }
			ss << amin << " #leq " << aleg << " < " << amax;
			legnames[Alpha] = ss.str();
			ss.str(string());
			ss << aprint << amin << "to" << amax;
			printnames[Alpha] = ss.str();
			printname += "_" + printnames[Alpha];

			ss.str(string());
			ss << ptmin << " #leq p_{T} < " << ptmax;
			legnames[Pt] = ss.str();
			ss.str(string());
			ss << "pt" << ptmin << "to" << ptmax;
			printnames[Pt] = ss.str();
			printname += "_" + printnames[Pt];
			
			ss.str(string());
			ss << etamin << " #leq |#eta| < " << etamax;
			legnames[Eta] = ss.str();
			ss.str(string());
			ss << "eta" << etamin << "to" << etamax;
			printnames[Eta] = ss.str();
			printname += "_" + printnames[Eta];
			
			legnames[Extra] = extra;
			printnames[Extra] = extra;
			
			//format response histogram
			hist->GetXaxis()->SetTitle("Asymmetry");
			hist->GetYaxis()->SetTitle("number of events");
			hist->SetLineWidth(2);
			hist->SetLineColor(color);
			
			//create log version of hist
			TH1F* loghist = (TH1F*)hist->Clone("log");
			for(int b = 0; b <= loghist->GetNbinsX()+1; ++b){
				if(loghist->GetBinContent(b)<=0) continue;
				loghist->SetBinError(b,loghist->GetBinError(b)/loghist->GetBinContent(b));
				loghist->SetBinContent(b,log(loghist->GetBinContent(b)));
			}
			
			//get mean value of alpha histo
			alpha_mean = h_alpha->GetMean();
			alpha_meanE = h_alpha->GetMeanError();
			
			//get values from histo
			mean  = loghist->GetMean();
			meanE = loghist->GetMeanError();
			rms   = loghist->GetRMS();
			rmsE  = loghist->GetRMSError();
			Nevents  = loghist->GetEntries();
			
			//fit histogram
			//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
			//ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
			//TVirtualFitter::SetMaxIterations(50000);
			TF1* logfit = new TF1("logtot",logcball,loghist->GetXaxis()->GetXmin(),loghist->GetXaxis()->GetXmax(),5);
			//assume tail starts at 2 sigma, peak is at 0
			logfit->SetParameters(hist->GetBinContent(1),0,rms,2,1.1);
			//limits on parameters: 0 < a < 10, 1 < n < 200
			logfit->SetParLimits(3,0,10);
			logfit->SetParLimits(4,1.01,200);
			loghist->Fit(logfit,"NQ");
			
			//get values from fit
			mu      = logfit->GetParameter(1);
			muE     = logfit->GetParError(1);
			sigma   = abs(logfit->GetParameter(2));
			sigmaE  = logfit->GetParError(2);
			a       = logfit->GetParameter(3);
			aE      = logfit->GetParError(3);
			n       = logfit->GetParameter(4);
			nE      = logfit->GetParError(4);
			//chi2ndf = logfit->GetChisquare()/logfit->GetNDF();
			
			//setup linear function
			gfit = new TF1("tot",cball,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),5);
			gfit->SetParameters(logfit->GetParameter(0),mu,sigma,a,n);
			//use linear function to compute chisquare with original histo
			//chi2ndf = hist->Chisquare(gfit)/logfit->GetNDF();
			
			//fit text
			fitnames.reserve(6);
			std::stringstream Nname, muname, sigmaname, aname, nname, chiname;
			//Nname << "N = " << Nevents; fitnames.push_back(Nname.str());
			muname.precision(3); muname << fixed << "#mu = " << mu << " #pm " << muE; fitnames.push_back(muname.str());
			sigmaname.precision(3); sigmaname << fixed << "#sigma = " << sigma << " #pm " << sigmaE; fitnames.push_back(sigmaname.str());
			aname.precision(3); aname << fixed << "a = " << a << " #pm " << aE; fitnames.push_back(aname.str());
			nname.precision(3); nname << fixed << "n = " << n << " #pm " << nE; fitnames.push_back(nname.str());
			chiname.precision(5); chiname << fixed << "#chi^{2}/ndf = " << chi2ndf; fitnames.push_back(chiname.str());
			
		}
	
		//accessors
		string GetDesc(qty q, bool print){ 
			if(print) return printnames[q];
			else return legnames[q];
		}
		//for plotting graphs of asym fit param vs. qty
		double GetX(qty q){
			switch(q){
				case Jet: return (double)jtype;
				case Alpha: return alpha_mean;
				case Pt: return pt;
				case Eta: return eta;
				default: return 0;
			}
		}
		double GetXerr(qty q){
			switch(q){
				case Alpha: return alpha_meanE;
				case Pt: return pt-ptmin;
				case Eta: return eta-etamin;
				default: return 0;
			}
		}
		double GetY(int p){
			switch(p){
				case 0: return mu;
				case 1: return sigma;
				case 2: return a;
				case 3: return n;
				default: return 0;
			}
		}
		double GetYerr(int p){
			switch(p){
				case 0: return muE;
				case 1: return sigmaE;
				case 2: return aE;
				case 3: return nE;
				default: return 0;
			}
		}
		
		//member variables
		alg jtype;
		alph atype;
		int year;
		double amin, amax, ptmin, ptmax, pt, etamin, etamax, eta, alpha_mean, alpha_meanE;
		string extra;
		Color_t color;
		vector<string> legnames, printnames, fitnames;
		string printname, cutname;
		TH1F* hist;
		TF1* gfit;
		double mean, meanE, rms, rmsE;
		double mu, muE, sigma, sigmaE, a, aE, n, nE, chi2ndf;
		int Nevents;
};

//----------------------------------------------------------------------------------------------------
//class to keep track of AsymFit results, what cuts they have in common, and to extrapolate fit params
//todo: make class that inherits from this one, but with vector asymextraps
//to plot extrap values vs. pt or eta
class KAsymExtrap {
	public:
		//constructor
		KAsymExtrap(Color_t color_=kBlack, int marker_=20) : color(color_), marker(marker_), asymfits(0), common(qtySize,true), q_varied(0) {
			for(int p = 0; p < 4; ++p){
				graph[p] = NULL;
				gfit[p] = NULL;
			}
		}
		
		//accessors
		void push_back(KAsymFit* s){
			asymfits.push_back(s);
			//keep track of common qtys among asymfits in this group
			if(asymfits.size()>1) compare(asymfits[0],s);
		}
		void compare(KAsymFit* s1, KAsymFit* s2){
			for(int q = 0; q < qtySize; q++){
				if(common[q]) common[q] = compare_qty(s1,s2,(qty)q);
			}
		}
		bool compare_qty(KAsymFit* s1, KAsymFit* s2, qty q){
			switch(q){
				case Jet: return s1->jtype==s2->jtype;
				case Alpha: return s1->alpha_mean==s2->alpha_mean;
				case Pt: return s1->ptmin==s2->ptmin && s1->ptmax==s2->ptmax;
				case Eta:  return s1->etamin==s2->etamin && s1->etamax==s2->etamax;
				case Extra: return s1->extra==s2->extra;
				default:   return true;
			}
		}
		//get descriptions from asymfits
		string GetDescAll(int s, bool print, bool incommon) {
			string desc("");
			for(int q = 0; q < qtySize; q++){
				string qdesc = GetDesc(s, (qty)q, print, incommon);
				if(qdesc.size()>0) {
					if(desc.size()>0 && print) desc += "_";
					else if(desc.size()>0) desc += ", ";
					desc += qdesc;
				}
			}
			return desc;
		}
		vector<string> GetDescList(int s, bool print, bool incommon){
			vector<string> descs;
			for(int q = 0; q < qtySize; q++){
				string qdesc = GetDesc(s, (qty)q, print, incommon);
				if(qdesc.size()>0) descs.push_back(qdesc);
			}
			return descs;
		}
		string GetDesc(int s, qty q, bool print, bool incommon) {
			if(asymfits.size()<s+1) return "";
			if(common[q]==incommon) return asymfits[s]->GetDesc(q, print);
			else return "";
		}
		//more clearly named instances from above
		string GetCommonDesc()  { return GetDescAll(0,false,true); }
		string GetCommonPrint() { return GetDescAll(0,true,true); }
		string GetVariedDesc(int s)  { return GetDescAll(s,false,false); }
		string GetVariedPrint(int s) { return GetDescAll(s,true,false); }
		vector<string> GetCommonDescList()  { return GetDescList(0,false,true); }
		vector<string> GetCommonPrintList() { return GetDescList(0,true,true); }
		vector<string> GetVariedDescList(int s)  { return GetDescList(s,false,false); }
		vector<string> GetVariedPrintList(int s) { return GetDescList(s,true,false); }
		//make param extrap graphs
		void MakeGraphs(bool fit){
			//check for q automatically
			q_varied = 0;
			for(int q = 0; q < qtySize; q++){
				if(!common[q]){
					if(q_varied==0) q_varied = q;
					else {
						cout << "Warning: multiple varied quantities in this group!" << endl;
					}
				}
			}
			if(q_varied==0) cout << "Warning: no varied quantities in this group!" << endl;
			
			for(int p = 0; p < 4; ++p){
				MakeGraph(p);
			}
		}
		void MakeGraph(int p){
			double* x = new double[asymfits.size()];
			double* xe = new double[asymfits.size()];
			double* y = new double[asymfits.size()];
			double* ye = new double[asymfits.size()];
			
			for(int s = 0; s < asymfits.size(); s++){
				x[s] = asymfits[s]->GetX((qty)q_varied);
				xe[s] = asymfits[s]->GetXerr((qty)q_varied);
				y[s] = asymfits[s]->GetY(p);
				ye[s] = asymfits[s]->GetYerr(p);
			}
			
			if(graph[p]) delete graph[p];
			graph[p] = new TGraphErrors(asymfits.size(),x,y,xe,ye);
			
			//axis titles
			switch((qty)q_varied){
				case Jet: graph[p]->GetXaxis()->SetTitle("Jet type"); break;
				case Alpha: 
					if(asymfits[0]->atype==Std) { graph[p]->GetXaxis()->SetTitle("#alpha"); }
					else if(asymfits[0]->atype==Par) { graph[p]->GetXaxis()->SetTitle("#alpha_{#parallel}"); }
					else if(asymfits[0]->atype==Perp) { graph[p]->GetXaxis()->SetTitle("#alpha_{#perp}"); }
					break;
				case Pt: graph[p]->GetXaxis()->SetTitle("p_{T} [GeV]"); break;
				case Eta:  graph[p]->GetXaxis()->SetTitle("#eta"); break;
				default:   graph[p]->GetXaxis()->SetTitle("");
			}
			switch(p){
				case 0: graph[p]->GetYaxis()->SetTitle("#mu");
				case 1: graph[p]->GetYaxis()->SetTitle("#sigma");
				case 2: graph[p]->GetYaxis()->SetTitle("a");
				case 3: graph[p]->GetYaxis()->SetTitle("n");
				default: graph[p]->GetYaxis()->SetTitle("");
			}
			
			//formatting
			graph[p]->SetMarkerStyle(marker);
			graph[p]->SetMarkerColor(color);
			graph[p]->SetLineColor(color);
			
			//linear fit
			if(gfit[p]) delete gfit[p];
			gfit[p] = new TF1("lin","pol1",graph[p]->GetXaxis()->GetXmin(),graph[p]->GetXaxis()->GetXmax());
			graph[p]->Fit(gfit[p],"NQ");
			
			//get values from fit
			p0[p]   = gfit[p]->GetParameter(0);
			p0E[p]  = gfit[p]->GetParError(0);
			chi2ndf[p] = gfit[p]->GetChisquare()/gfit[p]->GetNDF();
			
			//fit text
			fitnames[p].reserve(2);
			std::stringstream p0name, chiname;
			p0name.precision(3); p0name << fixed << graph[p]->GetYaxis()->GetTitle() << "_{0} = " << p0 << " #pm " << p0E; fitnames[p].push_back(p0name.str());
			chiname.precision(5); chiname << fixed << "#chi^{2}/ndf = " << chi2ndf[p]; fitnames[p].push_back(chiname.str());
		}
		
		//member variables
		Color_t color;
		int marker;
		vector<KAsymFit*> asymfits;
		vector<bool> common;
		TGraphErrors* graph[4];
		TF1* gfit[4];
		double p0[4], p0E[4], chi2ndf[4];
		int q_varied;
		vector<string> fitnames[4];
};

#endif
