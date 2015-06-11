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

#define nExtrapPars 5

using namespace std;
using namespace TMath;

enum qty { None=0, Jet=1, Alpha=2, Pt=3, Eta=4, Extra=5, qtySize=6 };
enum alg { Reco=0, Gen=1 };
enum alph { Std=0, Par=1, Perp=2 };

//------------------------------------------
//Symmetric Novosibirsk function
//following notation from CALICE in JINST 9 P01004
//parameters:
//N, mu, sigma, tau
//0,  1,     2,   3
Double_t novos(Double_t *x, Double_t *par){
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	//ensure sigma > 0
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	//ensure tau > 0 for right-side tail
	par[3] = Abs(par[3]);
	Double_t tau = par[3];
	Double_t arg = Abs(x[0]-mu)/sigma;
	
	//converges to gaussian as tau -> 0
	if(Abs(tau)< 1.e-7) return N*Exp(-Power(arg,2)/2);

	//sqrt(ln(4))
	static const Double_t xi = 1.177410022515475;
	Double_t lnarg = 1 + arg*SinH(xi*tau)/xi;
	
	//avoid small or negative arguments to log
	if(lnarg < 1.e-7) return 0.0;
	
	return N*Exp(-(Power(Log(lnarg)/tau,2)+Power(tau,2))/2);
}

//------------------------------------------
//Symmetric log(Novosibirsk) function
//following notation from CALICE in JINST 9 P01004
//parameters:
//N, mu, sigma, tau
//0,  1,     2,   3
Double_t lognovos(Double_t *x, Double_t *par){
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	//ensure sigma > 0
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	//ensure tau > 0 for right-side tail
	par[3] = Abs(par[3]);
	Double_t tau = par[3]; //tail param can be positive or negative?
	Double_t arg = Abs(x[0]-mu)/sigma;
	
	//converges to gaussian as tau -> 0
	if(Abs(tau) < 1.e-7) return Log(N) - Power(arg,2)/2;

	//sqrt(ln(4))
	static const Double_t xi = 1.177410022515475;
	Double_t lnarg = 1 + arg*SinH(xi*tau)/xi;
	
	//avoid small or negative arguments to log
	if(lnarg < 1.e-7) return 0.0;
	
	return Log(N) - (Power(Log(lnarg)/tau,2)+Power(tau,2))/2;
}

//-----------------------------------------------------------------------------------------
//class to keep track of set of cut values, corresponding histo, fit, various string labels
class KAsymFit {
	public:
		//constructor
		KAsymFit(TH1F* hist_, alg jtype_, alph atype_, double amin_, double amax_, double amean_, double ameanE_, double ptmin_, double ptmax_, double etamin_, double etamax_, string extra_="", Color_t color_=kBlack) :
			jtype(jtype_), atype(atype_), amin(amin_), amax(amax_), alpha_mean(amean_), alpha_meanE(ameanE_), ptmin(ptmin_), 
			ptmax(ptmax_), pt((ptmin+ptmax)/2), etamin(etamin_), etamax(etamax_), eta((etamin+etamax)/2), extra(extra_), color(color_),
			legnames(qtySize,""), printnames(qtySize,""), cutname(""), printname(""),
			hist(hist_), nfit(NULL), gfit(NULL), mean(0), meanE(0), rms(0), rmsE(0), Nevents(0),
			norm(0), normE(0), mu(0), muE(0), sigma(0), sigmaE(0), tau(0), tauE(0), chi2ndf(0),
			gnorm(0), gnormE(0), gmu(0), gmuE(0), gsigma(0), gsigmaE(0), gchi2ndf(0)
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
			
			//get values from histo
			mean  = loghist->GetMean();
			meanE = loghist->GetMeanError();
			rms   = loghist->GetRMS();
			rmsE  = loghist->GetRMSError();
			Nevents  = loghist->GetEntries();
			
			//fit histogram
			//use "smart fit" approach from JECs
			TF1* lognfit = new TF1("logtot",lognovos,loghist->GetXaxis()->GetXmin(),loghist->GetXaxis()->GetXmax(),4);
			//assume peak is at 0, arbitrary tail param
			lognfit->SetParameters(hist->GetBinContent(1),0,rms,0.5);
			//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
			//ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
			//TVirtualFitter::SetMaxIterations(50000);
			chi2ndf = -1;
			for(int f = 0; f < 30; ++f){
				loghist->Fit(lognfit,"NQ");
				
				double ctmp = lognfit->GetChisquare()/lognfit->GetNDF();
				if(ctmp < chi2ndf || chi2ndf == -1){
					//get values from fit
					norm    = lognfit->GetParameter(0);
					normE   = lognfit->GetParError(0);
					mu      = lognfit->GetParameter(1);
					muE     = lognfit->GetParError(1);
					sigma   = abs(lognfit->GetParameter(2));
					sigmaE  = lognfit->GetParError(2);
					tau     = abs(lognfit->GetParameter(3));
					tauE    = lognfit->GetParError(3);
					chi2ndf = lognfit->GetChisquare()/lognfit->GetNDF();
				}
			}
			
			//setup linear function
			nfit = new TF1("tot",novos,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),4);
			nfit->SetParameters(norm,mu,sigma,tau);
			//use linear function to compute chisquare with original histo
			//chi2ndf = hist->Chisquare(gfit)/logfit->GetNDF();
			
			//fit text
			fitnames.reserve(6);
			std::stringstream Nname, muname, sigmaname, tauname, chiname;
			//Nname << "N = " << Nevents; fitnames.push_back(Nname.str());
			muname.precision(3); muname << fixed << "#mu = " << mu << " #pm " << muE; fitnames.push_back(muname.str());
			sigmaname.precision(3); sigmaname << fixed << "#sigma = " << sigma << " #pm " << sigmaE; fitnames.push_back(sigmaname.str());
			tauname.precision(3); tauname << fixed << "#tau = " << tau << " #pm " << tauE; fitnames.push_back(tauname.str());
			chiname.precision(5); chiname << fixed << "#chi^{2}/ndf = " << chi2ndf; fitnames.push_back(chiname.str());
			
			//now do the "classic" gaussian core fit for comparison
			//iteration 1: fit to range determined by mean and RMS
			gfit = new TF1("gsn","gaus",max(0.,mean-2.5*rms),mean+2.5*rms);
			hist->Fit(gfit,"NQR");
			//iteration 2: fit to range determined by mu and sigma from iter 1
			gfit->SetRange(max(0.,gfit->GetParameter(1)-2*gfit->GetParameter(2)),gfit->GetParameter(1)+2*gfit->GetParameter(2));
			hist->Fit(gfit,"NQR");
			//get values from fit
			gnorm    = gfit->GetParameter(0);
			gnormE   = gfit->GetParError(0);
			gmu      = gfit->GetParameter(1);
			gmuE     = gfit->GetParError(1);
			gsigma   = abs(gfit->GetParameter(2));
			gsigmaE  = gfit->GetParError(2);
			gchi2ndf = gfit->GetChisquare()/gfit->GetNDF();
			//fit text
			std::stringstream gmuname, gsigmaname, gchiname;
			gmuname.precision(3); gmuname << fixed << "#mu = " << gmu << " #pm " << gmuE; fitnames.push_back(gmuname.str());
			gsigmaname.precision(3); gsigmaname << fixed << "#sigma = " << gsigma << " #pm " << gsigmaE; fitnames.push_back(gsigmaname.str());
			gchiname.precision(5); gchiname << fixed << "#chi^{2}/ndf = " << gchi2ndf; fitnames.push_back(gchiname.str());
			
			//todo: add tail measurements and uncertainties...
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
				case 2: return tau;
				case 3: return gmu;
				case 4: return gsigma;
				default: return 0;
			}
		}
		double GetYerr(int p){
			switch(p){
				case 0: return muE;
				case 1: return sigmaE;
				case 2: return tauE;
				case 3: return gmuE;
				case 4: return gsigmaE;
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
		TF1 *nfit, *gfit;
		double mean, meanE, rms, rmsE;
		double norm, normE, mu, muE, sigma, sigmaE, tau, tauE, chi2ndf;
		double gnorm, gnormE, gmu, gmuE, gsigma, gsigmaE, gchi2ndf;
		int Nevents;
};

//----------------------------------------------------------------------------------------------------
//class to keep track of AsymFit results, what cuts they have in common, and to extrapolate fit params
//todo: make class that inherits from this one, but with vector asymextraps
//to plot extrap values vs. pt or eta
class KAsymExtrap {
	public:
		//constructor
		KAsymExtrap(Color_t color_=kBlack, int marker_=20) : color(color_), marker(marker_), asymfits(0), asymextraps(0), common(qtySize,true), q_varied(0) {
			for(int p = 0; p < nExtrapPars; ++p){
				graph[p] = NULL;
				gfit[p] = NULL;
				ymax[p] = 0;
				ymin[p] = 1e10;
			}
		}
		
		//accessors
		virtual void push_back(KAsymFit* s){
			asymfits.push_back(s);
			//keep track of common qtys among asymfits in this group
			if(asymfits.size()>1) compare(asymfits[0],s);
		}
		virtual void compare(KAsymFit* s1, KAsymFit* s2){
			for(int q = 0; q < qtySize; q++){
				if(common[q]) common[q] = compare_qty(s1,s2,(qty)q);
			}
		}
		virtual bool compare_qty(KAsymFit* s1, KAsymFit* s2, qty q){
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
		virtual string GetDescAll(int s, bool print, bool incommon) {
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
		virtual vector<string> GetDescList(int s, bool print, bool incommon){
			vector<string> descs;
			for(int q = 0; q < qtySize; q++){
				string qdesc = GetDesc(s, (qty)q, print, incommon);
				if(qdesc.size()>0) descs.push_back(qdesc);
			}
			return descs;
		}
		virtual string GetDesc(int s, qty q, bool print, bool incommon) {
			if(asymfits.size()<s+1) return "";
			if(common[q]==incommon) return asymfits[s]->GetDesc(q, print);
			else return "";
		}
		//more clearly named instances from above
		virtual string GetCommonDesc()  { return GetDescAll(0,false,true); }
		virtual string GetCommonPrint() { return GetDescAll(0,true,true); }
		virtual string GetVariedDesc(int s)  { return GetDescAll(s,false,false); }
		virtual string GetVariedPrint(int s) { return GetDescAll(s,true,false); }
		virtual vector<string> GetCommonDescList()  { return GetDescList(0,false,true); }
		virtual vector<string> GetCommonPrintList() { return GetDescList(0,true,true); }
		virtual vector<string> GetVariedDescList(int s)  { return GetDescList(s,false,false); }
		virtual vector<string> GetVariedPrintList(int s) { return GetDescList(s,true,false); }
		//make param extrap graphs
		virtual void MakeGraphs(double xmin, double xmax){
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
			if(q_varied==0) { cout << "Warning: no varied quantities in this group!" << endl; }
			
			for(int p = 0; p < nExtrapPars; ++p){
				MakeGraph(p,xmin,xmax);
			}
		}
		virtual void MakeGraph(int p, double xmin, double xmax){
			double* x = new double[asymfits.size()];
			double* xe = new double[asymfits.size()];
			double* y = new double[asymfits.size()];
			double* ye = new double[asymfits.size()];
			
			for(int s = 0; s < asymfits.size(); s++){
				x[s] = asymfits[s]->GetX((qty)q_varied);
				xe[s] = asymfits[s]->GetXerr((qty)q_varied);
				y[s] = asymfits[s]->GetY(p);
				ye[s] = asymfits[s]->GetYerr(p);
				if(y[s]>ymax[p]) ymax[p] = y[s];
				if(y[s]<ymin[p]) ymin[p] = y[s];
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
				case 0: graph[p]->GetYaxis()->SetTitle("#mu"); break;
				case 1: graph[p]->GetYaxis()->SetTitle("#sigma"); break;
				case 2: graph[p]->GetYaxis()->SetTitle("#tau"); break;
				case 3: graph[p]->GetYaxis()->SetTitle("#mu_{gsn}"); break;
				case 4: graph[p]->GetYaxis()->SetTitle("#sigma_{gsn}"); break;
				default: graph[p]->GetYaxis()->SetTitle("");
			}
			
			//formatting
			graph[p]->SetMarkerStyle(marker);
			graph[p]->SetMarkerColor(color);
			graph[p]->SetLineColor(color);
			
			//linear fit
			if(gfit[p]) delete gfit[p];
			if(q_varied==Alpha) gfit[p] = new TF1("lin","pol1",xmin,xmax);
			else gfit[p] = new TF1("lin","pol1",graph[p]->GetXaxis()->GetXmin(),graph[p]->GetXaxis()->GetXmax());
			graph[p]->Fit(gfit[p],"NQ");
			
			//get values from fit
			p0[p]   = gfit[p]->GetParameter(0);
			p0E[p]  = gfit[p]->GetParError(0);
			chi2ndf[p] = gfit[p]->GetChisquare()/gfit[p]->GetNDF();
			if(p0[p]>ymax[p]) ymax[p] = p0[p];
			if(p0[p]<ymin[p]) ymin[p] = p0[p];
			
			//fit text
			fitnames[p].reserve(2);
			std::stringstream p0name, chiname;
			p0name.precision(3); p0name << fixed << graph[p]->GetYaxis()->GetTitle() << "_{0} = " << p0[p] << " #pm " << p0E[p]; fitnames[p].push_back(p0name.str());
			chiname.precision(5); chiname << fixed << "#chi^{2}/ndf = " << chi2ndf[p]; fitnames[p].push_back(chiname.str());
		}
		
		//member variables
		Color_t color;
		int marker;
		vector<KAsymFit*> asymfits;
		vector<KAsymExtrap*> asymextraps;
		vector<bool> common;
		TGraphErrors* graph[nExtrapPars];
		TF1* gfit[nExtrapPars];
		double p0[nExtrapPars], p0E[nExtrapPars], chi2ndf[nExtrapPars], ymax[nExtrapPars], ymin[nExtrapPars];
		int q_varied;
		vector<string> fitnames[nExtrapPars];
};

//----------------------------------------------------------------------------------------------------
//class to keep track of AsymExtrap results, what cuts they have in common, and look at trends
//inherits from KAsymExtrap
class KAsymTrend : public KAsymExtrap {
	public:
		//constructor
		KAsymTrend(Color_t color_=kBlack, int marker_=20) : KAsymExtrap(color_, marker_) { }
		
		//accessors
		void push_back(KAsymExtrap* s){
			asymextraps.push_back(s);
			//cout << s->GetCommonPrint() << ", " << s->GetVariedPrint(s->asymfits.size()-1) << endl;
			//keep track of common qtys among asymfits in this group
			if(asymextraps.size()>1) compare(asymextraps[0],s);
		}
		void compare(KAsymExtrap* s1, KAsymExtrap* s2){
			for(int q = 0; q < qtySize; q++){
				if(common[q]) {
					//if both extraps have the same varied quantity, it's "common" for the trend
					if(s1->q_varied==q){
						if(s2->q_varied==q) common[q] = true;
						else common[q] = false;
					}
					else {
						common[q] = compare_qty(s1->asymfits[0], s2->asymfits[0], (qty)q);
					}
				}
			}
		}
		//get descriptions from asymextraps
		string GetDesc(int s, qty q, bool print, bool incommon) {
			if(asymextraps.size()<s+1) return "";
			//second condition necessary to skip the var that's varied in each extrap
			if(common[q]==incommon){
				//special case
				if(q==asymextraps[s]->q_varied && q==Alpha){
					if(asymextraps[s]->asymfits[0]->atype==Std) return print ? "Alpha" : "#alpha";
					else if(asymextraps[s]->asymfits[0]->atype==Par) return print ? "Apar" : "#alpha_{#parallel}";
					else if(asymextraps[s]->asymfits[0]->atype==Perp) return print ? "Aperp" : "#alpha_{#perp}";
					else return "";
				}
				else return asymextraps[s]->GetDesc(0,q,print,incommon);
			}
			else return "";
		}
		//make param trend graphs
		void MakeGraph(int p, double xmin, double xmax){
			double* x = new double[asymextraps.size()];
			double* xe = new double[asymextraps.size()];
			double* y = new double[asymextraps.size()];
			double* ye = new double[asymextraps.size()];
			
			for(int s = 0; s < asymextraps.size(); s++){
				x[s] = asymextraps[s]->asymfits[0]->GetX((qty)q_varied);
				xe[s] = asymextraps[s]->asymfits[0]->GetXerr((qty)q_varied);
				y[s] = asymextraps[s]->p0[p];
				ye[s] = asymextraps[s]->p0E[p];
				if(y[s]>ymax[p]) ymax[p] = y[s];
				if(y[s]<ymin[p]) ymin[p] = y[s];
			}
			
			if(graph[p]) delete graph[p];
			graph[p] = new TGraphErrors(asymextraps.size(),x,y,xe,ye);
			
			//axis titles
			switch((qty)q_varied){
				case Jet: graph[p]->GetXaxis()->SetTitle("Jet type"); break;
				case Alpha: 
					if(asymextraps[0]->asymfits[0]->atype==Std) { graph[p]->GetXaxis()->SetTitle("#alpha"); }
					else if(asymextraps[0]->asymfits[0]->atype==Par) { graph[p]->GetXaxis()->SetTitle("#alpha_{#parallel}"); }
					else if(asymextraps[0]->asymfits[0]->atype==Perp) { graph[p]->GetXaxis()->SetTitle("#alpha_{#perp}"); }
					break;
				case Pt: graph[p]->GetXaxis()->SetTitle("p_{T} [GeV]"); break;
				case Eta:  graph[p]->GetXaxis()->SetTitle("#eta"); break;
				default:   graph[p]->GetXaxis()->SetTitle("");
			}
			switch(p){
				case 0: graph[p]->GetYaxis()->SetTitle("#mu"); break;
				case 1: graph[p]->GetYaxis()->SetTitle("#sigma"); break;
				case 2: graph[p]->GetYaxis()->SetTitle("#tau"); break;
				case 3: graph[p]->GetYaxis()->SetTitle("#mu_{gsn}"); break;
				case 4: graph[p]->GetYaxis()->SetTitle("#sigma_{gsn}"); break;
				default: graph[p]->GetYaxis()->SetTitle("");
			}
			
			//formatting
			graph[p]->SetMarkerStyle(marker);
			graph[p]->SetMarkerColor(color);
			graph[p]->SetLineColor(color);
		}
		
		//member variables
};

#endif