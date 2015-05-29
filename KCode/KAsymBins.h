//custom headers
#include "Analysis/KCode/KParser.h"

//STL headers
#include <vector>
#include <string>
#include <fstream>

using namespace std;

//define multidimensional vector types
typedef vector<string>  svec1;
typedef vector<svec1>   svec2;
typedef vector<svec2>   svec3;
typedef vector<svec3>   svec4;
typedef vector<svec4>   svec5;
typedef vector<TH1F*>  hvec1;
typedef vector<hvec1>   hvec2;
typedef vector<hvec2>   hvec3;
typedef vector<hvec3>   hvec4;
typedef vector<hvec4>   hvec5;

class KAsymBins {
	public:
		//constructor
		KAsymBins(string input="input/input_bins.txt") : parsed(false) {
			//must always have an option map
			globalOpt = new OptionMap();
			
			//parse input
			string intype;
			string line;
			ifstream infile(input.c_str());
			if(infile.is_open()){
				while(getline(infile,line)){
					//skip commented lines
					if(line[0]=='#') continue;
					//skip blank lines
					if(line.size()<2) continue;
					
					//check for carriage returns (not allowed)
					if(line[line.size()-1]=='\r') {
						cout << "Carriage return detected. Please run:" << endl;
						cout << "dos2unix " << input << endl;
						cout << "and then try again." << endl;
						return;
					}
					
					//check for input type
					if(line.compare(0,6,"OPTION")==0) { intype = "OPTION"; continue; }
					
					//otherwise, process line according to input type
					if(intype=="OPTION") KParser::processOption(line,globalOpt);
				}
			}
			parsed = true;
			
			//get vectors
			globalOpt->Get("pt_bins",pt);
			globalOpt->Get("eta_bins",eta);
			globalOpt->Get("alpha_bins",alpha);
			globalOpt->Get("alpha_types",atype);
			globalOpt->Get("jet_types",jtype);
		}
		
		//helper functions
		void MakeHistos(bool only_names=false){
			//initialize name vectors
			asym_incl_names = asym_excl_names = alpha_incl_names = alpha_excl_names = svec5(jtype.size(),svec4(atype.size(),svec3(pt.size(),svec2(eta.size(),svec1(alpha.size(),"")))));
			
			//initialize histo vectors if requested
			if(!only_names){
				asym_incl = asym_excl = alpha_incl = alpha_excl = hvec5(jtype.size(),hvec4(atype.size(),hvec3(pt.size(),hvec2(eta.size(),hvec1(alpha.size(),NULL)))));
			}
			
			//axes: jtype, atype, pt, eta, alpha
			for(int ijt = 0; ijt < jtype.size(); ++ijt){
				for(int iat = 0; iat < atype.size(); ++iat){
					for(int ipt = 0; ipt < pt.size()-1; ++ipt){
						for(int iet = 0; iet < eta.size()-1; ++iet){
							for(int ial = 0; ial < alpha.size()-1; ++ial){
								//make names in easily parsable format
								stringstream nss;
								nss << "jtype" << ijt << "_" << "atype" << iat << "_" << "pt" << ipt << "_" << "eta" << iet << "_" << "alpha" << ial;
								string ntmp = nss.str();
								asym_incl_names[ijt][iat][ipt][iet][ial] = "asym__" + ntmp + "in";
								asym_excl_names[ijt][iat][ipt][iet][ial] = "asym__" + ntmp + "ex";
								alpha_incl_names[ijt][iat][ipt][iet][ial] = "alpha__" + ntmp + "in";
								alpha_excl_names[ijt][iat][ipt][iet][ial] = "alpha__" + ntmp + "ex";
								
								//make histos if requested
								if(!only_names){
									asym_incl[ijt][iat][ipt][iet][ial] = new TH1F(asym_incl_names[ijt][iat][ipt][iet][ial].c_str(),"",50,0.0,1.0);
									asym_excl[ijt][iat][ipt][iet][ial] = new TH1F(asym_excl_names[ijt][iat][ipt][iet][ial].c_str(),"",50,0.0,1.0);
									alpha_incl[ijt][iat][ipt][iet][ial] = new TH1F(alpha_incl_names[ijt][iat][ipt][iet][ial].c_str(),"",50,alpha[ial],alpha.back());
									alpha_excl[ijt][iat][ipt][iet][ial] = new TH1F(alpha_excl_names[ijt][iat][ipt][iet][ial].c_str(),"",50,alpha[ial],alpha[ial+1]);
								}
							}
						}
					}
				}
			}
		}
		
		//member variables
		bool parsed;
		OptionMap* globalOpt;
		vector<double> pt, eta, alpha;
		vector<int> atype, jtype;
		svec5 asym_incl_names, asym_excl_names, alpha_incl_names, alpha_excl_names;
		hvec5 asym_incl, asym_excl, alpha_incl, alpha_excl;
};