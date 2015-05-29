#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "KCode/TreeClass.C"

using namespace std;

//to run interactively:
//root -b -q -l KHistDriver.C+
void KHistDriver(){
	TFile* file = TFile::Open("tree_dijet/tree_QCD.root");
	TTree* tree = (TTree*)file->Get("tree");
	TreeClass t(tree);
	t.Loop();
}