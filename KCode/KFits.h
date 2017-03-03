#ifndef KFITS_H
#define KFITS_H

//STL headers
#include <iostream>

//ROOT headers
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TLine.h>

//custom headers
#include "KFit.h"
#include "Analysis/KCode/KMath.h"

using namespace std;
using namespace TMath;

//----------------------------------------------
//fit that selects the better of pol1 or pol0
class KLinOrConstFit : public KFitGraph {
	public:
		//constructors
		KLinOrConstFit() : KFitGraph() {}
		KLinOrConstFit(string name_, OptionMap* localOpt_) : KFitGraph(name_,localOpt_) {}
		//destructor
		virtual ~KLinOrConstFit() {}
		
		//functions
		virtual void Fit() {
			if(!dist) {
				cout << "Error: fit " << name << " has no dist to fit!" << endl;
				return;
			}
			
			//better of pol0 vs pol1 fit
			TGraph* graph = dist->GetDist();
			TF1* gfit0 = new TF1("const","pol0",graph->GetXaxis()->GetXmin(),graph->GetXaxis()->GetXmax());
			graph->Fit(gfit0,"NQ");
			TF1* gfit1 = new TF1("lin","pol1",graph->GetXaxis()->GetXmin(),graph->GetXaxis()->GetXmax());
			graph->Fit(gfit1,"NQ");
			
			//choose better fit (or fit specified by options)
			if(localOpt->Get("const",false)) { fit = gfit0; delete gfit1; }
			else if(localOpt->Get("lin",false)) { fit = gfit1; delete gfit0; }
			else if(gfit0->GetChisquare()/gfit0->GetNDF() < gfit1->GetChisquare()/gfit1->GetNDF()) { fit = gfit0; delete gfit1; }
			else { fit = gfit1; delete gfit0; }
			
			//make fit quantities
			vector<string> dname_fields;
			KParser::process(dist->GetSourceQty(),':',dname_fields);
			string dpname = dname_fields.size()==3 ? dname_fields[2] : dist->GetName();
			KFitQuantity* q_p0 = new KFitQuantity(dpname+"0",dist->LegName()+"_{0}",dist->GetName()+"0",dist->AxisName()+"_{0}",fit->GetParameter(0),fit->GetParError(0));
			StoreFitQty(q_p0);
			KFitQuantity* q_chi2ndf = new KFitQuantity("chi2ndf","#chi^{2}/ndf","chi2ndf","#chi^{2}/ndf",fit->GetChisquare()/fit->GetNDF());
			StoreFitQty(q_chi2ndf);
		}
};

//----------------------------------------------
//"classic" R+S gaussian fit
class KGaussianCoreFit : public KFitHisto {
	public:
		//constructors
		KGaussianCoreFit() : KFitHisto() {}
		KGaussianCoreFit(string name_, OptionMap* localOpt_) : KFitHisto(name_,localOpt_) {}
		//destructor
		virtual ~KGaussianCoreFit() {}
		
		//integrals and derivatives
		double pdf(double x, double* p){
			return p[0]*Exp(-0.5*Power((x-p[1])/p[2],2));
		}
		double cdf(double x, double* p){
			return p[0]*p[2]*Sqrt(Pi()/2)*Erf((x-p[1])/(Sqrt(2)*p[2]));
		}
		double dcdf_dN(double x, double* p){
			return cdf(x,p)/p[0];
		}
		double dcdf_dmu(double x, double* p){
			return -pdf(x,p);
		}
		double dcdf_dsigma(double x, double* p){
			return cdf(x,p)/p[2] - pdf(x,p)*(x-p[1])/p[2];
		}
		double dcdf_dsigma_v2(double x, double* p){ //when x = a*sigma
			return cdf(x,p)/p[2] + pdf(x,p)*p[1]/p[2];
		}
		
		//functions
		virtual void Fit() {
			if(!dist) {
				cout << "Error: fit " << name << " has no dist to fit!" << endl;
				return;
			}
			
			//get histos
			TH1F* hist = dist->GetDist();
			vector<TH1F*> hsplit;
			dist->GetLocalOpt()->Get("hsplit",hsplit);
			vector<float> weights;
			dist->GetLocalOpt()->Get("weights",weights);
			
			//get values from histo
			double mean  = hist->GetMean();
			double meanE = hist->GetMeanError();
			double rms   = hist->GetRMS();
			double rmsE  = hist->GetRMSError();
			int Nevents  = hist->GetEntries();
			
			//now do the "classic" gaussian core fit for comparison
			//iteration 1: fit to range determined by mean and RMS
			fit = new TF1("gsn","gaus",max(0.,mean-2.5*rms),mean+2.5*rms);
			hist->Fit(fit,"NQR");
			//iteration 2: fit to range determined by mu and sigma from iter 1
			double sigma_iter1 = fit->GetParameter(2);
			fit->SetRange(max(0.,fit->GetParameter(1)-2*fit->GetParameter(2)),fit->GetParameter(1)+2*fit->GetParameter(2));
			TFitResultPtr result = hist->Fit(fit,"NQRS");
			
			//make fit quantities
			KFitQuantity* q_mu = new KFitQuantity("mu","#mu","gmu","#mu",fit->GetParameter(1),fit->GetParError(1));
			StoreFitQty(q_mu);
			KFitQuantity* q_sigma = new KFitQuantity("sigma","#sigma","gsigma","#sigma",fit->GetParameter(2),fit->GetParError(2));
			StoreFitQty(q_sigma);
			KFitQuantity* q_chi2ndf = new KFitQuantity("chi2ndf","#chi^{2}/ndf","gchi2ndf","#chi^{2}/ndf",fit->GetChisquare()/fit->GetNDF());
			StoreFitQty(q_chi2ndf);
			
			//option to skip tail frac calc
			if(!localOpt->Get("do_ftail",true)) return;
			
			//todo: split tail fraction MC measurement and data measurement into separate functions
			//idea: if hsplit not available, assume data
			//then dependency check for other fit result (store result in localopt manually in macro)
			//if check succeeds, calc data/mc SF
			
			//calculate Gaussian part of integral
			double* pars = fit->GetParameters();
			double* errs = fit->GetParErrors();
			TMatrixDSym cov = result->GetCovarianceMatrix();
			double Ibound_b = 1.0;
			double Ibound_a = 2.5*q_sigma->GetVal();
			double Ibound_a_low = Ibound_a - 2.5*q_sigma->GetErrLow();
			double Ibound_a_up = Ibound_a + 2.5*q_sigma->GetErrHigh();
			double Igaus_b = cdf(Ibound_b,pars);
			double Igaus_a = cdf(Ibound_a,pars);
			double Igaus = Igaus_b - Igaus_a;
			
			//propagate uncertainties in Gaussian part of integral
			double dIgaus_dN = dcdf_dN(Ibound_b,pars) - dcdf_dN(Ibound_a,pars);
			double dIgaus_dmu = dcdf_dmu(Ibound_b,pars) - dcdf_dmu(Ibound_a,pars);
			double dIgaus_dsigma = dcdf_dsigma(Ibound_b,pars) - dcdf_dsigma_v2(Ibound_a,pars); //account for extra dependence on sigma in bound
			double var_Igaus = Power(dIgaus_dN*errs[0],2) + Power(dIgaus_dmu*errs[1],2) + Power(dIgaus_dsigma*errs[2],2);
			//add in covariance if available
			if(cov.GetNrows()>=3 && cov.GetNcols()>=3) var_Igaus += 2*( dIgaus_dN*dIgaus_dmu*cov[0][1] + dIgaus_dN*dIgaus_dsigma*cov[0][2] + dIgaus_dmu*dIgaus_dsigma*cov[1][2] );
			double std_Igaus = Sqrt(var_Igaus);
			
			//denominator for Gaussian tail fraction: integral of histo
			//split into [0,2*sigma_iter1], [2*sigma_iter1,1] : error on first part is 100% correlated with Gaussian fit param errors, second part is 100% uncorrelated
			double std_Ihist_low = 0;
			int sigma_iter1_bin = hist->FindBin(2*sigma_iter1);
			hist->IntegralAndError(0,sigma_iter1_bin,std_Ihist_low);
			double std_Ihist = 0;
			double Ihist = hist->IntegralAndError(0,hist->GetNbinsX()+1,std_Ihist);
			double var_Ihist = Power(std_Ihist,2);
			double ftail_gaus = Igaus/Ihist;
			double dftail_gaus_dgaus = 1./Ihist;
			double dftail_gaus_dhist = -Igaus/Power(Ihist,2);
			double var_ftail_gaus = var_Igaus*Power(dftail_gaus_dgaus,2) + var_Ihist*Power(dftail_gaus_dhist,2) + 2*dftail_gaus_dgaus*dftail_gaus_dhist*std_Igaus*std_Ihist_low;
			//todo: put Ibound_a, ftail_gaus, var_ftail_gaus into DIST localOpt so the novos fit can also access them

			//set boundary of integral
			//options: gsigma, gsigma_0 (extrapolated to alpha=0), gsigma_C (0 <= alpha < 0.05)
			//currently using gsigma
			//todo: test other options
			//propagate uncertainty on bound as independent error on fraction, sigma_p/N
			int sigma_bin = hist->FindBin(Ibound_a);
			int sigma_bin_low = hist->FindBin(Ibound_a_low);
			int sigma_bin_up = hist->FindBin(Ibound_a_up);
			
			//tail fraction measurement from MC (use event weights for errors)
			//for each pT^hat bin with different xsec weight:
			//calculate the error terms by hand, then sum up
			double fmc_numer, fmc_numer_low, fmc_numer_up, fmc_denom;
			fmc_numer = fmc_numer_low = fmc_numer_up = fmc_denom = 0;
			double p_err_low, p_err_up, f_err_low, f_err_up;
			p_err_low = p_err_up = f_err_low = f_err_up = 0.;
			for(int iwe = 0; iwe < weights.size(); ++iwe){
				TH1F* htmp = hsplit[iwe];
				if(!htmp || htmp->GetEntries()==0) continue;
				
				//get total, pass, fail
				int N_i = htmp->Integral(0,htmp->GetNbinsX()+1);
				int p_i = htmp->Integral(sigma_bin,htmp->GetNbinsX()+1);
				int p_i_low = htmp->Integral(sigma_bin_low,htmp->GetNbinsX()+1);
				int p_i_up = htmp->Integral(sigma_bin_up,htmp->GetNbinsX()+1);
				int f_i = N_i - p_i;

				//add to efficiency factors
				fmc_numer += weights[iwe]*p_i;
				fmc_numer_low += weights[iwe]*p_i_low;
				fmc_numer_up += weights[iwe]*p_i_up;
				fmc_denom += weights[iwe]*N_i;
				
				//calculate errors
				p_err_low += Power(weights[iwe],2)*Power(KMath::PoissonErrorLow(p_i),2);
				p_err_up += Power(weights[iwe],2)*Power(KMath::PoissonErrorUp(p_i),2);
				f_err_low += Power(weights[iwe],2)*Power(KMath::PoissonErrorLow(f_i),2);
				f_err_up += Power(weights[iwe],2)*Power(KMath::PoissonErrorUp(f_i),2);
			}
			//finish error calculation
			//NB: this is a lazy way of propagating asymmetric errors
			double fmc = fmc_numer/fmc_denom;
			double fmcElow = Sqrt( Power(1-fmc,2)*p_err_low + Power(fmc,2)*f_err_low + (localOpt->Get("propagate_sigma_bin",true) ? Power(fmc_numer-fmc_numer_low,2) : 0) )/fmc_denom;
			double fmcEup = Sqrt( Power(1-fmc,2)*p_err_up + Power(fmc,2)*f_err_up + (localOpt->Get("propagate_sigma_bin",true) ? Power(fmc_numer-fmc_numer_up,2) : 0) )/fmc_denom;
			//subtract ftail_gaus (neglecting cov(fmc,ftail_gaus) uncertainty term)
			if(localOpt->Get("subtract_gaus",true)){
				fmc = fmc - ftail_gaus;
				fmcElow = Sqrt(Power(fmcElow,2)+var_ftail_gaus);
				fmcEup = Sqrt(Power(fmcEup,2)+var_ftail_gaus);
			}
			//make fit quantity
			KFitQuantity* q_fmc = new KFitQuantity("fmc","#it{f}_{mc}","gfmc","#it{f}_{mc}",fmc,fmcElow,fmcEup);
			StoreFitQty(q_fmc);
			
			//tail fraction measurement from "data" : no bin weights, simple efficiency
			//todo: for tail fraction from novos fit, clone histo with Poisson bin errors
			//todo: account for prescale factor based on pTave bin?
			double fdata = 0;
			double fdataElow = 0;
			double fdataEup = 0;
			if(hist->GetEntries()>0){
				double N_data = round(hist->Integral(0,hist->GetNbinsX()+1));
				double p_data = round(hist->Integral(sigma_bin,hist->GetNbinsX()+1));
				double p_data_low = round(hist->Integral(sigma_bin_low,hist->GetNbinsX()+1));
				double p_data_up = round(hist->Integral(sigma_bin_up,hist->GetNbinsX()+1));
				double f_data = N_data - p_data;
				fdata = p_data/N_data;
				fdataElow = Sqrt( Power(1-fdata,2)*Power(KMath::PoissonErrorLow((int)p_data),2) + Power(fdata,2)*Power(KMath::PoissonErrorLow((int)f_data),2) + (localOpt->Get("propagate_sigma_bin",true) ? Power(p_data-p_data_low,2) : 0) )/N_data;
				fdataEup = Sqrt( Power(1-fdata,2)*Power(KMath::PoissonErrorUp((int)p_data),2) + Power(fdata,2)*Power(KMath::PoissonErrorUp((int)f_data),2) + (localOpt->Get("propagate_sigma_bin",true) ? Power(p_data-p_data_up,2) : 0) )/N_data;
				//subtract ftail_gaus (neglecting cov(fdata,ftail_gaus) uncertainty term)
				if(localOpt->Get("subtract_gaus",true)){
					fdata = fdata - ftail_gaus;
					fdataElow = Sqrt(Power(fdataElow,2)+var_ftail_gaus);
					fdataEup = Sqrt(Power(fdataEup,2)+var_ftail_gaus);
				}
			}
			//make fit quantity
			KFitQuantity* q_fdata = new KFitQuantity("fdata","#it{f}_{data}","gfdata","#it{f}_{data}",fdata,fdataElow,fdataEup);
			StoreFitQty(q_fdata);
				
			//scale factor
			double sf = 0;
			double sfElow = 0;
			double sfEup = 0;
			if(fmc!=0){
				sf = fdata/fmc;
				sfElow = sf*Sqrt(Power(fdataElow,2)/Power(fdata,2) + Power(fmcElow,2)/Power(fmc,2));
				sfEup = sf*Sqrt(Power(fdataEup,2)/Power(fdata,2) + Power(fmcEup,2)/Power(fmc,2));
			}
			else { //nan protection
				sf = 1.0;
				sfElow = 0.0;
				sfEup = 0.0;
			}
			//make fit quantity
			KFitQuantity* q_sf = new KFitQuantity("sf","#it{s}_{data/mc}","gsf","#it{s}_{data/mc}",sf,sfElow,sfEup);
			StoreFitQty(q_sf);
		}
		
		virtual void Draw(TPad* pad) {
			//plot the classic gaussian fit, extended to the full range
			fit->SetRange(dist->GetDist()->GetXaxis()->GetXmin(),dist->GetDist()->GetXaxis()->GetXmax());
			
			KFit::Draw(pad);
			
			//plot the tail start line
			Double_t bnd = fit->GetParameter(1) + fit->GetParameter(2)*2.5;
			double aymin = pad->GetLogy() ? Power(10,pad->GetUymin()) : pad->GetUymin();
			double aymax = pad->GetLogy() ? Power(10,pad->GetUymax()) : pad->GetUymax();
			TLine* aline = new TLine(bnd,aymin,bnd,aymax);
			aline->SetLineColor(color);
			aline->SetLineStyle(line);
			aline->SetLineWidth(width);
			aline->Draw("same");
		}
};

//----------------------------------------------
//right-sided crystal ball fit
class KCrystalBallFit : public KFitHisto {
	public:
		//constructors
		KCrystalBallFit() : KFitHisto() {}
		KCrystalBallFit(string name_, OptionMap* localOpt_) : KFitHisto(name_,localOpt_) {}
		//destructor
		virtual ~KCrystalBallFit() {}
		
		//fit function
		//parameters:
		//N, mu, sigma, a, n
		//0,  1,     2, 3, 4
		double cball(double *x, double *par){
			//ensure sigma > 0 and a > 0
			double N = par[0]; //let N float
			double mu = par[1];
			par[2] = Abs(par[2]);
			double sigma = par[2];
			double a = par[3];
			par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
			double n = par[4];
			double arg = (x[0]-mu)/sigma;
			
			//left tail
			//right tail
			if(arg >= a){
				return N*Power(n/a,n)*Exp(-Power(a,2)/2)*Power(n/a-a+arg,-n);
			}
			//core
			else{
				return N*Exp(-Power(arg,2)/2);
			}
		}
		//log(fit function)
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
		
		//integrals and derivatives
		double cdfgsn(double x, double* p){
			return p[0]*p[2]*Sqrt(Pi()/2)*Erf((x-p[1])/(Sqrt(2)*p[2]));
		}
		double cdfpwr(double x, double*p){
			return 1/(1-p[4])*Power(p[4]/p[3],p[4])*Exp(-0.5*Power(p[3],2))*Power(p[4]/p[3]-p[3]+(x-p[1])/p[2],-p[4]+1);
		}
		
		//functions
		virtual void Fit() {
			if(!dist) {
				cout << "Error: fit " << name << " has no dist to fit!" << endl;
				return;
			}
			
			//get histos
			TH1F* hist = dist->GetDist();
			
			//get values from histo
			double mean  = hist->GetMean();
			double meanE = hist->GetMeanError();
			double rms   = hist->GetRMS();
			double rmsE  = hist->GetRMSError();
			int Nevents  = hist->GetEntries();
			
			//crystal ball fit from member function
			fit = new TF1("cball",this,&KCrystalBallFit::cball,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),5);
			
			//option for log fit
			if(localOpt->Get("logfit",true)){
				//create log version of hist
				TH1F* loghist = (TH1F*)hist->Clone("log");
				for(int b = 0; b <= loghist->GetNbinsX()+1; ++b){
					if(loghist->GetBinContent(b)<=0) continue;
					loghist->SetBinError(b,loghist->GetBinError(b)/loghist->GetBinContent(b));
					loghist->SetBinContent(b,log(loghist->GetBinContent(b)));
				}
				
				//get values from loghist
				double logmean  = loghist->GetMean();
				double logmeanE = loghist->GetMeanError();
				double logrms   = loghist->GetRMS();
				double logrmsE  = loghist->GetRMSError();
				
				TF1* logfit = new TF1("logcball",this,&KCrystalBallFit::logcball,loghist->GetXaxis()->GetXmin(),loghist->GetXaxis()->GetXmax(),5);
				//initial values
				logfit->SetParameters(hist->GetBinContent(1),0,logrms,2.5,1.1);
				//fix the tail start to 2.5 sigma
				//todo: make configurable
				logfit->FixParameter(3,2.5);
				logfit->SetParLimits(4,1.01,200);
				//use smart fit procedure
				TFitResultPtr result = SmartFit(loghist,logfit);
				
				//put values in linear fit
				for(int p = 0; p < fit->GetNpar(); ++p){
					fit->SetParameter(p,logfit->GetParameter(p));
					fit->SetParError(p,logfit->GetParError(p));
				}
				fit->SetChisquare(logfit->GetChisquare());
				fit->SetNDF(logfit->GetNDF());
			}
			else {
				//initial values
				fit->SetParameters(hist->GetBinContent(1),0,rms,2.5,1.1);
				//fix the tail start to 2.5 sigma
				//todo: make configurable
				fit->FixParameter(3,2.5);
				//use smart fit procedure
				TFitResultPtr result = SmartFit(hist,fit);
			}
			
			//make fit quantities
			KFitQuantity* q_mu = new KFitQuantity("mu","#mu","cmu","#mu",fit->GetParameter(1),fit->GetParError(1));
			StoreFitQty(q_mu);
			KFitQuantity* q_sigma = new KFitQuantity("sigma","#sigma","csigma","#sigma",fit->GetParameter(2),fit->GetParError(2));
			StoreFitQty(q_sigma);
			KFitQuantity* q_n = new KFitQuantity("n","n","cn","n",fit->GetParameter(4),fit->GetParError(4));
			StoreFitQty(q_n);
			KFitQuantity* q_chi2ndf = new KFitQuantity("chi2ndf","#chi^{2}/ndf","cchi2ndf","#chi^{2}/ndf",fit->GetChisquare()/fit->GetNDF());
			StoreFitQty(q_chi2ndf);
			
			//todo: tail fraction with uncertainty propagation
		}
		
		virtual void Draw(TPad* pad) {
			KFit::Draw(pad);
			
			//plot the tail start line
			Double_t bnd = fit->GetParameter(1) + fit->GetParameter(2)*fit->GetParameter(3);
			double aymin = pad->GetLogy() ? Power(10,pad->GetUymin()) : pad->GetUymin();
			double aymax = pad->GetLogy() ? Power(10,pad->GetUymax()) : pad->GetUymax();
			TLine* aline = new TLine(bnd,aymin,bnd,aymax);
			aline->SetLineColor(color);
			aline->SetLineStyle(line+1);
			aline->SetLineWidth(width);
			aline->Draw("same");
			
			//plot the gaussian part
			TF1* gsn = new TF1("gsn","gaus",bnd,dist->GetDist()->GetXaxis()->GetXmax());
			gsn->SetParameters(fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));
			gsn->SetLineColor(color);
			gsn->SetLineWidth(line+1);
			gsn->SetLineStyle(width);
			gsn->Draw("same");		
		}
};

//addition to KParser to create fits
namespace KParser {
	KFit* processFit(KNamed* tmp){
		KFit* ftmp = 0;
		string name = tmp->first;
		OptionMap* omap = tmp->second;
		
		//construct if known fit
		if(name=="linorconst") ftmp = new KLinOrConstFit(name,omap);
		else if(name=="gaussiancore") ftmp = new KGaussianCoreFit(name,omap);
		else if(name=="crystalball") ftmp = new KCrystalBallFit(name,omap);
		
		if(!ftmp) cout << "Input error: unknown fit " << name << " will be skipped." << endl;
		
		return ftmp;
	}
}


#endif