//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 28 14:55:58 2015 by ROOT version 5.34/18
// from TTree tree/selected observables, dijet
// found on file: tree_dijet/tree_QCD.root
//////////////////////////////////////////////////////////

#ifndef TreeClass_h
#define TreeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Float_t         weight;
   vector<TLorentzVector> *RecoJet;
   Double_t        RecoAlpha;
   Double_t        RecoApar;
   Double_t        RecoAperp;
   Double_t        RecoDijetPt;
   Double_t        RecoDijetDeltaPhi;
   Double_t        RecoAsymmetry;
   vector<TLorentzVector> *GenJet;
   Double_t        GenAlpha;
   Double_t        GenApar;
   Double_t        GenAperp;
   Double_t        GenDijetPt;
   Double_t        GenDijetDeltaPhi;
   Double_t        GenAsymmetry;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_RecoJet;   //!
   TBranch        *b_recoalpha;   //!
   TBranch        *b_recoapar;   //!
   TBranch        *b_recoaperp;   //!
   TBranch        *b_recodijetpt;   //!
   TBranch        *b_recodijetdphi;   //!
   TBranch        *b_recoasymmetry;   //!
   TBranch        *b_GenJet;   //!
   TBranch        *b_genalpha;   //!
   TBranch        *b_genapar;   //!
   TBranch        *b_genaperp;   //!
   TBranch        *b_gendijetpt;   //!
   TBranch        *b_gendijetdphi;   //!
   TBranch        *b_genasymmetry;   //!

   TreeClass(TTree *tree=0);
   virtual ~TreeClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeClass_cxx
TreeClass::TreeClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tree_dijet/tree_QCD.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tree_dijet/tree_QCD.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

TreeClass::~TreeClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   RecoJet = 0;
   GenJet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("RecoJet", &RecoJet, &b_RecoJet);
   fChain->SetBranchAddress("RecoAlpha", &RecoAlpha, &b_recoalpha);
   fChain->SetBranchAddress("RecoApar", &RecoApar, &b_recoapar);
   fChain->SetBranchAddress("RecoAperp", &RecoAperp, &b_recoaperp);
   fChain->SetBranchAddress("RecoDijetPt", &RecoDijetPt, &b_recodijetpt);
   fChain->SetBranchAddress("RecoDijetDeltaPhi", &RecoDijetDeltaPhi, &b_recodijetdphi);
   fChain->SetBranchAddress("RecoAsymmetry", &RecoAsymmetry, &b_recoasymmetry);
   fChain->SetBranchAddress("GenJet", &GenJet, &b_GenJet);
   fChain->SetBranchAddress("GenAlpha", &GenAlpha, &b_genalpha);
   fChain->SetBranchAddress("GenApar", &GenApar, &b_genapar);
   fChain->SetBranchAddress("GenAperp", &GenAperp, &b_genaperp);
   fChain->SetBranchAddress("GenDijetPt", &GenDijetPt, &b_gendijetpt);
   fChain->SetBranchAddress("GenDijetDeltaPhi", &GenDijetDeltaPhi, &b_gendijetdphi);
   fChain->SetBranchAddress("GenAsymmetry", &GenAsymmetry, &b_genasymmetry);
   Notify();
}

Bool_t TreeClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeClass_cxx
