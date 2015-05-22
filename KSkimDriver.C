#include <TROOT.h>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "KCode/KSkimManager.h"

using namespace std;

//to recompile:
//root -b -q -l 'KSkimDriver.C++()'
//to run interactively:
//root -b -q -l 'KSkimDriver.C+("input/input_selection.txt","QCD","dijet"," root://cmseos.fnal.gov//store/user/pedrok/SUSY2015/Phys14_QCD_Pt-binned_PUPPI","tree")'
void KSkimDriver(string input="", string setname="", string selTypes="", string indir="", string outdir="tree"){
	if(input.size()==0){
		cout << "Recompiled KSkimDriver, exiting." << endl;
		return;
	}
	KSkimManager k(input,setname,selTypes,indir,outdir);
	k.Skim();
}