#ifndef KTREND_H
#define KTREND_H

//STL headers
#include <string>
#include <sstream>

//ROOT headers
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>

//custom headers
#include "Analysis/KCode/KMap.h"
#include "Analysis/KCode/KPlot.h"
#include "Analysis/KCode/KLegend.h"
#include "KQuantities.h"
#include "KDist.h"
#include "KFits.h"

using namespace std;
using namespace TMath;

class KTrend {
	public:
		//constructor
		KTrend() :
			localOpt(new OptionMap()), globalOpt(new OptionMap()), name(""), color(kBlack), marker(20), q_varied(-1), fit_quantities(new FitQtyMap()) {}
		KTrend(OptionMap* localOpt_, OptionMap* globalOpt_) : 
			localOpt(localOpt_ ? localOpt_ : new OptionMap()), globalOpt(globalOpt_ ? globalOpt_ : new OptionMap()), name(""), color(kBlack), marker(20), q_varied(-1), fit_quantities(new FitQtyMap()) {}
		//destructor
		~KTrend() {}
		
		//initialization
		void Initialize(){
			//step 0: get basic attributes
			localOpt->Get("name",name);
			//bool setColor = localOpt->Get("color",color);
			//bool setMarker = localOpt->Get("marker",marker);
			
			//step 1: check for children and quantities
			bool hasChildren = localOpt->Get("children",children);
			vector<string> quantity_names;
			bool hasQuantities = localOpt->Get("quantities",quantity_names);
			
			//step 2a: if children, check for common/uncommon/varied qtys
			if(hasChildren && children.size()>0){
				//get quantities from first child
				quantities = children[0]->quantities;
				common.resize(quantities.size(),1);
				
				CompareChildren();
			}
			//step 2b: if no children, initialize qtys
			else if(hasQuantities && quantity_names.size()>0){
				MakeQuantities(quantity_names);
				common.resize(quantities.size(),1);
			}
			//otherwise, no info for initialization, quit
			else return;
			
			//step 3: get/create dists (which create and execute their own fits)
			vector<KNamed*> distLines;
			bool hasDists = localOpt->Get("distLines",distLines);
			if(hasDists && distLines.size()>0) MakeDists(distLines);
			else cout << "Warning: no dists specified for trend " << name << "!" << endl;
		}
		void CompareChildren(){
			//possible values: 0 = varied (!common), 1 = common, 2 = distributed (varied among all children -> common for this trend)
			
			//loop over children, compare vs child 0
			for(unsigned c = 0; c < children.size(); ++c){
				if(c==0) continue;
				for(unsigned q = 0; q < quantities.size(); ++q){
					//so far common, needs to be checked
					if(common[q]) {
						//check for distributed status:
						//a qty which is varied for one child but not another is automatically varied for the trend
						if(children[0]->q_varied==q) common[q] = children[c]->q_varied==q ? 2 : 0;
						//also check for distributed status further down the line
						else if(children[0]->common[q]==2) common[q] = children[c]->common[q]==2 ? 2 : 0;
						//otherwise, compare quantities directly
						else common[q] = (*(children[0]->quantities[q]) == *(children[c]->quantities[q]));
					}
				}
			}
			
			//find the varied quantity
			for(unsigned q = 0; q < quantities.size(); ++q){
				if(!common[q]){
					if(q_varied==-1) q_varied = q;
					else {
						cout << "Warning: multiple varied quantities in trend " << name << "!" << endl;
					}
				}
			}
			if(q_varied==-1) { cout << "Warning: no varied quantities in trend " << name << "!" << endl; }
		}
		void MakeQuantities(vector<string>& quantity_names){
			quantities.reserve(quantity_names.size());
			for(unsigned q = 0; q < quantity_names.size(); ++q){
				KQuantity* qtmp = KParser::processQuantity(new KNamed(quantity_names[q],localOpt));
				if(qtmp) quantities.push_back(qtmp);
			}
		}
		void MakeDists(vector<KNamed*>& distLines){
			dists.reserve(distLines.size());
			for(unsigned d = 0; d < distLines.size(); ++d){
				KDist* dtmp = KParser::processDist(distLines[d]);
				if(!dtmp) continue;
				
				//check if the dist needs to be built
				
				//simplest case: existing object, nothing to do
				if(dtmp->GetSource()==KDist::Object) {}
				//construct using TTree::Draw
				else if(dtmp->GetSource()==KDist::Tree){
					//todo: different options for Histo or Graph (from V1, V2)
					if(dtmp->GetType()==KDist::Graph) {
						
					}
					else if(dtmp->GetType()==KDist::Histo) {
						//optionally, use KPlot::CreateHist() to make histo from parameters?
					}
				}
				//construct using KFitQuantities from children
				else if(dtmp->GetSource()==KDist::Children){
					//only allow Graph
					if(dtmp->GetType()==KDist::Graph) {
						//check availability of fit quantity
						KFitQuantity* fq0 = children[0]->fit_quantities->Get(dtmp->GetSourceQty());
						if(fq0){
							//set options for dist
							dtmp->GetLocalOpt()->Set("legname",fq0->LegName());
							dtmp->GetLocalOpt()->Set("axisname",fq0->AxisName());
							dtmp->GetLocalOpt()->Set("xtitle",children[0]->quantities[q_varied]->AxisName());
							dtmp->GetLocalOpt()->Set("ytitle",fq0->AxisName());
							
							//graph arrays
							double* x = new double[children.size()];
							double* xel = new double[children.size()];
							double* xeu = new double[children.size()];
							double* y = new double[children.size()];
							double* yel = new double[children.size()];
							double* yeu = new double[children.size()];
							
							//loop over children
							//todo: allow x,y to come from both quantities and fit quantities?
							for(unsigned c = 0; c < children.size(); ++c){
								//x from quantity
								KQuantity* q = children[c]->quantities[q_varied];
								x[c] = q->GetVal();
								xel[c] = q->GetErrLow();
								xeu[c] = q->GetErrHigh();
								
								//y from fit quantity
								KFitQuantity* fq = children[c]->fit_quantities->Get(dtmp->GetSourceQty());
								y[c] = fq->GetVal();
								yel[c] = fq->GetErrLow();
								yeu[c] = fq->GetErrHigh();
							}
							
							//check for y-axis range if needed
							if(!(dtmp->GetLocalOpt()->Has("yfixed"))){
								double ymin = MinElement(children.size(),y);
								double ymax = MaxElement(children.size(),y);
								double yrange = ymax - ymin;
								ymin = ymin - yrange*0.1;
								ymax = ymax + yrange*0.1;
								dtmp->GetLocalOpt()->Set<double>("ymin",ymin);
								dtmp->GetLocalOpt()->Set<double>("ymax",ymax);
							}
							
							//make graph
							TGraphAsymmErrors* gtmp = new TGraphAsymmErrors(children.size(),x,y,xel,xeu,yel,yeu);
							//set in dist
							dtmp->GetLocalOpt()->Set("dist",gtmp);
							
							//assemble dist print name as xxx_vs_(q_varied)
							//using just qty part of dist name
							vector<string> dname_fields;
							KParser::process(dtmp->GetSourceQty(),':',dname_fields);
							string dpname = dname_fields.size()==3 ? dname_fields[2] : dtmp->GetName();
							string pname = dpname + "_vs_" + children[0]->quantities[q_varied]->PrintVsName();
							dtmp->GetLocalOpt()->Set("printname",pname);
						}
						else {
							cout << "Input error: cannot find fit quantity " << dtmp->GetSourceQty() << " in Children. This dist will not be initialized." << endl;
						}
					}
					else if(dtmp->GetType()==KDist::Histo) {
						cout << "Input error: cannot generate Histo from Children for dist " << dtmp->GetName() << ". This dist will not be initialized." << endl;
					}
				}
				//todo: add "construct as ratio of two distributions" case?
				//unknown source, skip
				else {
					cout << "Input error: unknown source " << dtmp->GetSource() << " for dist " << dtmp->GetName() << ". This dist will not be initialized." << endl;
				}
				
				//pass the color and marker information to the dist (if not set separately)
				//if(!dtmp->GetLocalOpt()->Has("color")) dtmp->GetLocalOpt()->Set("color",color);
				//if(!dtmp->GetLocalOpt()->Has("marker")) dtmp->GetLocalOpt()->Set("marker",marker);

				//everything is ready for initialization
				dtmp->Initialize();
				if(dtmp->Initialized()) {
					//get fit qtys from dist fits
					FitQtyMap* fqm = dtmp->GetFitQtys();
					fit_quantities->GetTable().insert(fqm->GetTable().begin(),fqm->GetTable().end());

					//store dist
					dists.push_back(dtmp);
				}
				else {
					cout << "Input error: dist " << dtmp->GetName() << " could not be initialized." << endl;
					delete dtmp;
				}
			}
		}

		//get descriptions
		string GetDesc(unsigned q, bool print, bool vs){
			if(q>=quantities.size()) return "";
			if(vs){
				if(print) return quantities[q]->PrintVsName();
				else return quantities[q]->LegVsName();
			}
			else{
				if(print) return quantities[q]->PrintName();
				else return quantities[q]->LegName();
			}
		}
		string GetDesc(int c, unsigned q, bool print, bool incommon) {
			if(c>=(int)children.size() || q>=quantities.size()) return "";
			//cast common to bool
			if((common[q]!=0)==incommon) {
				if(c<0) return GetDesc(q, print, common[q]==2);
				else return children[c]->GetDesc(q, print, common[q]==2);
			}
			else return "";
		}
		string GetDescAll(int c, bool print, bool incommon) {
			string desc("");
			for(unsigned q = 0; q < quantities.size(); ++q){
				string qdesc = GetDesc(c, q, print, incommon);
				if(qdesc.size()>0) {
					if(desc.size()>0){
						if(print) desc += "_";
						else desc += ", ";
					}
					desc += qdesc;
				}
			}
			return desc;
		}
		vector<string> GetDescList(int c, bool print, bool incommon){
			vector<string> descs;
			for(unsigned q = 0; q < quantities.size(); ++q){
				string qdesc = GetDesc(c, q, print, incommon);
				if(qdesc.size()>0) descs.push_back(qdesc);
			}
			return descs;
		}
		//more clearly named instances from above
		string GetCommonLeg()  { return GetDescAll(-1,false,true); }
		string GetCommonPrint() { return GetDescAll(-1,true,true); }
		string GetVariedLeg(int c)  { return GetDescAll(c,false,false); }
		string GetVariedPrint(int c) { return GetDescAll(c,true,false); }
		vector<string> GetCommonLegList()  { return GetDescList(-1,false,true); }
		vector<string> GetCommonPrintList() { return GetDescList(-1,true,true); }
		vector<string> GetVariedLegList(int c)  { return GetDescList(c,false,false); }
		vector<string> GetVariedPrintList(int c) { return GetDescList(c,true,false); }
		
		//drawing
		void Draw(bool doSelf=false, bool doSelfFits=false, bool doChildren=false, bool doChildrenFits=false, bool doPrint=false){
			if(!doSelf && !doChildren) return;
			
			if(doSelf){
				DrawLoop(true,doSelfFits,dists,doPrint);
			}
			
			if(doChildren && children.size()>0){
				DrawLoop(false,doChildrenFits,children[0]->dists,doPrint);
			}
		}
		void DrawLoop(bool doSelf, bool doFits, vector<KDist*>& thedists, bool doPrint){
			//loop over dists
			for(unsigned d = 0; d < thedists.size(); ++d){
				//get base histo from dist
				TH1F* hbase = thedists[d]->GetBaseHist();
				
				//check y-axis range
				if(doSelf){
					double ymin, ymax;
					if(thedists[d]->GetLocalOpt()->Get("ymin",ymin) && thedists[d]->GetLocalOpt()->Get("ymax",ymax)){
						hbase->GetYaxis()->SetRangeUser(ymin,ymax);
					}
				}
				else {
					double *ymins = new double[children.size()];
					double *ymaxs = new double[children.size()];
					for(unsigned c = 0; c < children.size(); ++c){
						children[c]->dists[d]->GetLocalOpt()->Get("ymin",ymins[c]);
						children[c]->dists[d]->GetLocalOpt()->Get("ymax",ymaxs[c]);
					}
					double ymin = MinElement(children.size(),ymins);
					double ymax = MaxElement(children.size(),ymaxs);
					hbase->GetYaxis()->SetRangeUser(ymin,ymax);
				}
				
				//make canvas/print name
				string oname = thedists[d]->PrintName();
				oname += "__" + GetCommonPrint();
				if(!doSelf){
					for(unsigned c = 0; c < children.size(); ++c){
						string varied = GetVariedPrint(c);
						if(varied.size()>0) oname += "__" + varied;
					}
				}

				//get plot options from dist
				OptionMap *plotLocalOpt, *plotGlobalOpt;
				if(!(thedists[d]->GetLocalOpt()->Get("plotLocalOpt",plotLocalOpt))) plotLocalOpt = new OptionMap();
				if(!(thedists[d]->GetLocalOpt()->Get("plotGlobalOpt",plotGlobalOpt))) plotGlobalOpt = new OptionMap();
				
				//get preamble text
				vector<string> preamble = GetCommonLegList();
				plotGlobalOpt->Set("extra_text",preamble);
				//preamble in panel 0 by default - ignore if panel balancing is enabled
				if(!plotGlobalOpt->Get("balance_panels",true)){
					vector<int> extra_text_panels(preamble.size(),0);
					plotGlobalOpt->Set("extra_text_panels",extra_text_panels);
				}
				
				//make plot (with options from dist)
				KPlot* plot = new KPlot(oname,plotLocalOpt,plotGlobalOpt);
				plot->Initialize(hbase);
				KLegend* kleg = plot->GetLegend();
				TCanvas* can = plot->GetCanvas();
				TPad* pad1 = plot->GetPad1();
				pad1->cd();
				
				//make sure legend has a reference histo (for tick sizes)
				if(thedists[d]->GetType()==KDist::Graph) kleg->AddHist(hbase);
				
				//legend entries
				if(doSelf){
					AddToLegend(doSelf,doFits,thedists[d],kleg);
				}
				else{
					//add children to legend
					for(unsigned c = 0; c < children.size(); ++c){
						AddToLegend(doSelf,doFits,children[c]->dists[d],kleg,c);
					}
				}
				
				//build legend
				KLegend::Horz horz = KLegend::hdefault;
				KLegend::Vert vert = KLegend::vdefault;
				thedists[d]->GetLocalOpt()->Get("leg_horz",horz);
				thedists[d]->GetLocalOpt()->Get("leg_vert",vert);
				if(vert==KLegend::vdefault) kleg->Build(horz);
				else kleg->Build(horz,vert);

				//draw blank histo for axes
				plot->DrawHist();
				
				//draw legend and text
				plot->DrawText();
				
				//draw dists & fits
				if(doSelf){
					DrawDist(doFits,thedists[d],pad1);
				}
				else{
					for(unsigned c = 0; c < children.size(); ++c){
						DrawDist(doFits,children[c]->dists[d],pad1);
					}
				}
				
				//finish drawing
				plot->GetHisto()->Draw("sameaxis"); //draw again so axes on top

				//print if required
				if(doPrint) PrintCanvas(oname,can);
			}
		}
		void AddToLegend(bool doSelf, bool doFits, KDist* thedist, KLegend* kleg, int c=-1){
			//make sure legend knows about this histo
			if(thedist->GetType()==KDist::Histo) kleg->AddHist(static_cast<KDistHisto*>(thedist)->GetDist());
			
			//add child dist to legend
			if(!doSelf){
				thedist->AddToLegend(kleg,GetVariedLeg(c));
			}
			
			if(doFits){
				for(unsigned f = 0; f < thedist->GetFits().size(); ++f){
					//second param enables "addToPrev" in KLegend
					//(needed to keep fits associated with dists for children in balance_panels case)
					thedist->GetFits()[f]->AddToLegend(kleg,!doSelf);
				}
			}
		}
		void DrawDist(bool doFits, KDist* thedist, TPad* pad){
			thedist->Draw(pad);
			if(doFits){
				for(unsigned f = 0; f < thedist->GetFits().size(); ++f){
					thedist->GetFits()[f]->Draw(pad);
				}
			}
		}
		void PrintCanvas(string oname, TCanvas* can){
			//get print params from global options
			string psuff = "png";
			globalOpt->Get("psuff",psuff);
			string pdir = ".";
			globalOpt->Get("pdir",pdir);
			can->Print((pdir+"/"+oname+"."+psuff).c_str(),psuff.c_str());
		}
		
		//accessors
		OptionMap* GetLocalOpt() { return localOpt; }
		OptionMap* GetGlobalOpt() { return globalOpt; }
		
	protected:
		//member variables
		OptionMap* localOpt;
		OptionMap* globalOpt;
		string name;
		Color_t color;
		int marker;
		vector<KTrend*> children;
		vector<int> common;
		int q_varied;
		vector<KQuantity*> quantities;
		vector<KDist*> dists;
		FitQtyMap* fit_quantities;
		
};

#endif