//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 18 10:19:08 2014 by ROOT version 5.32/00
// from TTree pfTree/dijet tree
// found on file: /mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root
//////////////////////////////////////////////////////////

#ifndef pfCand_h
#define pfCand_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class pfCand {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nPFpart;
   Int_t           pfId[7082];   //[nPFpart]
   Float_t         pfPt[7082];   //[nPFpart]
   Float_t         pfVsPt[7082];   //[nPFpart]
   Float_t         pfVsPtInitial[7082];   //[nPFpart]
   Float_t         pfArea[7082];   //[nPFpart]
   Float_t         pfEta[7082];   //[nPFpart]
   Float_t         pfPhi[7082];   //[nPFpart]
   Float_t         vn[5][15];
   Float_t         psin[5][15];
   Float_t         sumpt[15];

   // List of branches
   TBranch        *b_nPFpart;   //!
   TBranch        *b_pfId;   //!
   TBranch        *b_pfPt;   //!
   TBranch        *b_pfVsPt;   //!
   TBranch        *b_pfVsPtInitial;   //!
   TBranch        *b_pfArea;   //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_vn;   //!
   TBranch        *b_vpsi;   //!
   TBranch        *b_sumpt;   //!

   pfCand(TString infile, TTree *tree=0);
   virtual ~pfCand();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Close();

   TFile * f;
};

#endif

#ifdef pfCand_cxx
pfCand::pfCand(TString infile, TTree *tree) 
{
   f = TFile::Open(infile);
   tree = (TTree*) f->Get("pfcandAnalyzer/pfTree");

   Init(tree);
}

pfCand::~pfCand()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pfCand::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pfCand::LoadTree(Long64_t entry)
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

void pfCand::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
   fChain->SetBranchAddress("pfId", pfId, &b_pfId);
   fChain->SetBranchAddress("pfPt", pfPt, &b_pfPt);
   fChain->SetBranchAddress("pfVsPt", pfVsPt, &b_pfVsPt);
   fChain->SetBranchAddress("pfVsPtInitial", pfVsPtInitial, &b_pfVsPtInitial);
   fChain->SetBranchAddress("pfArea", pfArea, &b_pfArea);
   fChain->SetBranchAddress("pfEta", pfEta, &b_pfEta);
   fChain->SetBranchAddress("pfPhi", pfPhi, &b_pfPhi);
   fChain->SetBranchAddress("vn", vn, &b_vn);
   fChain->SetBranchAddress("psin", psin, &b_vpsi);
   fChain->SetBranchAddress("sumpt", sumpt, &b_sumpt);
   Notify();
}

Bool_t pfCand::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pfCand::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pfCand::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void pfCand::Close()
{
  f->Close();
}
#endif // #ifdef pfCand_cxx
