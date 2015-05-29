//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 22 11:46:34 2015 by ROOT version 5.34/18
// from TTree t/t
// found on file: root://cmseos.fnal.gov//store/user/pedrok/SUSY2015/Phys14_QCD_Pt-binned_PUPPI/Pt-30to50.root
//////////////////////////////////////////////////////////

#ifndef NtupleClass_h
#define NtupleClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<int>     *npus;
   vector<float>   *tnpus;
   vector<int>     *bxns;
   Float_t         rho;
   Float_t         beta;
   Float_t         betaStar;
   Long64_t        npv;
   Long64_t        run;
   Long64_t        lumi;
   Long64_t        evt;
   Int_t           nref;
   vector<int>     *refrank;
   vector<int>     *refpdgid;
   vector<int>     *refpdgid_algorithmicDef;
   vector<int>     *refpdgid_physicsDef;
   vector<float>   *refe;
   vector<float>   *refpt;
   vector<float>   *refeta;
   vector<float>   *refphi;
   vector<float>   *refm;
   vector<float>   *refy;
   vector<float>   *refdrjt;
   vector<float>   *refarea;
   vector<float>   *jte;
   vector<float>   *jtpt;
   vector<float>   *jteta;
   vector<float>   *jtphi;
   vector<float>   *jtm;
   vector<float>   *jty;
   vector<float>   *jtjec;
   vector<float>   *jtarea;
   vector<float>   *jtchf;
   vector<float>   *jtnhf;
   vector<float>   *jtnef;
   vector<float>   *jtcef;
   vector<float>   *jtmuf;
   vector<float>   *jthfhf;
   vector<float>   *jthfef;
   UChar_t         nmu;
   Float_t         weight;

   // List of branches
   TBranch        *b_npus;   //!
   TBranch        *b_tnpus;   //!
   TBranch        *b_bxns;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_betaStar;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_refrank;   //!
   TBranch        *b_refpdgid;   //!
   TBranch        *b_refpdgid_algorithmicDef;   //!
   TBranch        *b_refpdgid_physicsDef;   //!
   TBranch        *b_refe;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refm;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_jte;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtjec;   //!
   TBranch        *b_jtarea;   //!
   TBranch        *b_jtchf;   //!
   TBranch        *b_jtnhf;   //!
   TBranch        *b_jtnef;   //!
   TBranch        *b_jtcef;   //!
   TBranch        *b_jtmuf;   //!
   TBranch        *b_jthfhf;   //!
   TBranch        *b_jthfef;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_weight;   //!

   NtupleClass(TTree *tree=0);
   virtual ~NtupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleClass_cxx
NtupleClass::NtupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/pedrok/SUSY2015/Phys14_QCD_Pt-binned_PUPPI/Pt-30to50.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmseos.fnal.gov//store/user/pedrok/SUSY2015/Phys14_QCD_Pt-binned_PUPPI/Pt-30to50.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

NtupleClass::~NtupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleClass::LoadTree(Long64_t entry)
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

void NtupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   npus = 0;
   tnpus = 0;
   bxns = 0;
   refrank = 0;
   refpdgid = 0;
   refpdgid_algorithmicDef = 0;
   refpdgid_physicsDef = 0;
   refe = 0;
   refpt = 0;
   refeta = 0;
   refphi = 0;
   refm = 0;
   refy = 0;
   refdrjt = 0;
   refarea = 0;
   jte = 0;
   jtpt = 0;
   jteta = 0;
   jtphi = 0;
   jtm = 0;
   jty = 0;
   jtjec = 0;
   jtarea = 0;
   jtchf = 0;
   jtnhf = 0;
   jtnef = 0;
   jtcef = 0;
   jtmuf = 0;
   jthfhf = 0;
   jthfef = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npus", &npus, &b_npus);
   fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
   fChain->SetBranchAddress("bxns", &bxns, &b_bxns);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("betaStar", &betaStar, &b_betaStar);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("refrank", &refrank, &b_refrank);
   fChain->SetBranchAddress("refpdgid", &refpdgid, &b_refpdgid);
   fChain->SetBranchAddress("refpdgid_algorithmicDef", &refpdgid_algorithmicDef, &b_refpdgid_algorithmicDef);
   fChain->SetBranchAddress("refpdgid_physicsDef", &refpdgid_physicsDef, &b_refpdgid_physicsDef);
   fChain->SetBranchAddress("refe", &refe, &b_refe);
   fChain->SetBranchAddress("refpt", &refpt, &b_refpt);
   fChain->SetBranchAddress("refeta", &refeta, &b_refeta);
   fChain->SetBranchAddress("refphi", &refphi, &b_refphi);
   fChain->SetBranchAddress("refm", &refm, &b_refm);
   fChain->SetBranchAddress("refy", &refy, &b_refy);
   fChain->SetBranchAddress("refdrjt", &refdrjt, &b_refdrjt);
   fChain->SetBranchAddress("refarea", &refarea, &b_refarea);
   fChain->SetBranchAddress("jte", &jte, &b_jte);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtm", &jtm, &b_jtm);
   fChain->SetBranchAddress("jty", &jty, &b_jty);
   fChain->SetBranchAddress("jtjec", &jtjec, &b_jtjec);
   fChain->SetBranchAddress("jtarea", &jtarea, &b_jtarea);
   fChain->SetBranchAddress("jtchf", &jtchf, &b_jtchf);
   fChain->SetBranchAddress("jtnhf", &jtnhf, &b_jtnhf);
   fChain->SetBranchAddress("jtnef", &jtnef, &b_jtnef);
   fChain->SetBranchAddress("jtcef", &jtcef, &b_jtcef);
   fChain->SetBranchAddress("jtmuf", &jtmuf, &b_jtmuf);
   fChain->SetBranchAddress("jthfhf", &jthfhf, &b_jthfhf);
   fChain->SetBranchAddress("jthfef", &jthfef, &b_jthfef);
   fChain->SetBranchAddress("nmu", &nmu, &b_nmu);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t NtupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupleClass_cxx
