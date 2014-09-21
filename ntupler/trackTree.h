//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 27 10:28:20 2013 by ROOT version 5.32/00
// from TTree trackTree/v1
// found on file: /mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet50_HydjetDrum_v27_mergedV1.root
//////////////////////////////////////////////////////////

#ifndef trackTree_h
#define trackTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "HltTree.C"
#include "HiTree.C"
#include "t.C"
#include "skimTree.C"
#include "genPart.C"
#include "hi.C"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class trackTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   // skimTree		  *fskim;
   // HltTree		  *fhlt;
   // t 	  *fjet;
   // HiTree		  *fhi;
   // genPart   *fgen;
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nEv;
   //Int_t           nLumi;
   //Int_t           nBX;
   //Int_t           nRun;
   //Int_t           nv;
   //Float_t         vx[2];   //[nv]
   //Float_t         vy[2];   //[nv]
   //Float_t         vz[2];   //[nv]
   Int_t           nTrk;
   Float_t         trkPt[10000];   //[nTrk]
   Float_t         trkPtError[10000];   //[nTrk]
   Int_t           trkNHit[10000];   //[nTrk]
   //Int_t           trkNlayer[10000];   //[nTrk]
   Float_t         trkEta[10000];   //[nTrk]
   Float_t         trkPhi[10000];   //[nTrk]
   //Int_t           trkCharge[10000];   //[nTrk]
   Bool_t          highPurity[10000];   //[nTrk]
   //Bool_t          highPuritySetWithPV[10000];   //[nTrk]
   //Float_t         trkChi2[10000];   //[nTrk]
   //Float_t         trkNdof[10000];   //[nTrk]
   Float_t         trkDxy1[10000];   //[nTrk]
   Float_t         trkDxyError1[10000];   //[nTrk]
   Float_t         trkDz1[10000];   //[nTrk]
   Float_t         trkDzError1[10000];   //[nTrk]
   Bool_t          trkFake[10000];   //[nTrk]
   Float_t         trkAlgo[10000];   //[nTrk]
   //Int_t           pfType[10000];   //[nTrk]
   //Float_t         pfCandPt[10000];   //[nTrk]
   //Float_t         pfSumEcal[10000];   //[nTrk]
   //Float_t         pfSumHcal[10000];   //[nTrk]
   Float_t         trkStatus[10000];   //[nTrk]
   Int_t           nParticle;
   //Float_t         pStatus[20000];   //[nParticle]
   //Float_t         pPId[20000];   //[nParticle]
   Float_t         pEta[20000];   //[nParticle]
   Float_t         pPhi[20000];   //[nParticle]
   Float_t         pPt[20000];   //[nParticle]
   //Float_t         pAcc[20000];   //[nParticle]
   //Float_t         pAccPair[20000];   //[nParticle]
   //Float_t         pNRec[20000];   //[nParticle]
   //Int_t           pNHit[20000];   //[nParticle]
   Float_t         mtrkPt[20000];   //[nParticle]
   Float_t         mtrkPtError[20000];   //[nParticle]
   //Int_t           mtrkNHit[20000];   //[nParticle]
   //Int_t           mtrkNlayer[20000];   //[nParticle]
   //Int_t           mtrkNlayer3D[20000];   //[nParticle]
   Int_t           mtrkQual[20000];   //[nParticle]
   //Float_t         mtrkChi2[20000];   //[nParticle]
   //Float_t         mtrkNdof[20000];   //[nParticle]
   Float_t         mtrkDz1[20000];   //[nParticle]
   Float_t         mtrkDzError1[20000];   //[nParticle]
   Float_t         mtrkDxy1[20000];   //[nParticle]
   Float_t         mtrkDxyError1[20000];   //[nParticle]
   //Float_t         mtrkAlgo[20000];   //[nParticle]
   //Int_t           mtrkPfType[20000];   //[nParticle]
   //Float_t         mtrkPfCandPt[20000];   //[nParticle]
   //Float_t         mtrkPfSumEcal[20000];   //[nParticle]
   //Float_t         mtrkPfSumHcal[20000];   //[nParticle]

   // List of branches
   TBranch        *b_nEv;   //!
   //TBranch        *b_nLumi;   //!
   //TBranch        *b_nBX;   //!
   //TBranch        *b_nRun;   //!
   //TBranch        *b_nv;   //!
   //TBranch        *b_vx;   //!
   //TBranch        *b_vy;   //!
   //TBranch        *b_vz;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   //TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   //TBranch        *b_trkCharge;   //!
   TBranch        *b_highPurity;   //!
   //TBranch        *b_highPuritySetWithPV;   //!
   //TBranch        *b_trkChi2;   //!
   //TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   //TBranch        *b_pfType;   //!
   //TBranch        *b_pfCandPt;   //!
   //TBranch        *b_pfSumEcal;   //!
   //TBranch        *b_pfSumHcal;   //!
   TBranch        *b_trkStatus;   //!
   TBranch        *b_nParticle;   //!
   //TBranch        *b_pStatus;   //!
   //TBranch        *b_pPId;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
   //TBranch        *b_pAcc;   //!
   //TBranch        *b_pAccPair;   //!
   //TBranch        *b_pNRec;   //!
   //TBranch        *b_pNHit;   //!
   TBranch        *b_mtrkPt;   //!
   TBranch        *b_mtrkPtError;   //!
   //TBranch        *b_mtrkNHit;   //!
   //TBranch        *b_mtrkNlayer;   //!
   //TBranch        *b_mtrkNlayer3D;   //!
   TBranch        *b_mtrkQual;   //!
   //TBranch        *b_mtrkChi2;   //!
   //TBranch        *b_mtrkNdof;   //!
   TBranch        *b_mtrkDz1;   //!
   TBranch        *b_mtrkDzError1;   //!
   TBranch        *b_mtrkDxy1;   //!
   TBranch        *b_mtrkDxyError1;   //!
   //TBranch        *b_mtrkAlgo;   //!
   //TBranch        *b_mtrkPfType;   //!
   //TBranch        *b_mtrkPfCandPt;   //!
   //TBranch        *b_mtrkPfSumEcal;   //!
   //TBranch        *b_mtrkPfSumHcal;   //!

   trackTree(TString infile,TTree *tree=0);

   virtual ~trackTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual Int_t    GetEntriesFast();
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void             Close();

   TFile *f;
};

#endif

#ifdef trackTree_cxx
trackTree::trackTree(TString infile,TTree *tree)
{
   f = TFile::Open(infile);
   std::cout<<infile.Data()<<std::endl;
   tree = (TTree*) f->Get("anaTrack/trackTree");
   // fhi = new HiTree(infile);
   // fhlt = new HltTree(infile);
   // fskim = new skimTree(infile);
   // fjet = new t(infile);
   // fgen = new genPart(infile);
   Init(tree);
}

trackTree::~trackTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trackTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
   // return fhi->GetEntry(entry);
   // return fhlt->GetEntry(entry);
   // return fskim->GetEntry(entry);
   // return fjet->GetEntry(entry);
   // return fgen->GetEntry(entry);
}

Int_t trackTree::GetEntriesFast()
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntriesFast();
}

Long64_t trackTree::LoadTree(Long64_t entry)
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

void trackTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   //fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   //fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   //fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   //fChain->SetBranchAddress("nv", &nv, &b_nv);
   //fChain->SetBranchAddress("vx", vx, &b_vx);
   //fChain->SetBranchAddress("vy", vy, &b_vy);
   //fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   //fChain->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   //fChain->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   //fChain->SetBranchAddress("highPuritySetWithPV", highPuritySetWithPV, &b_highPuritySetWithPV);
   //fChain->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   //fChain->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   fChain->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   fChain->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   fChain->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   fChain->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   fChain->SetBranchAddress("trkFake", trkFake, &b_trkFake);
   fChain->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
   //fChain->SetBranchAddress("pfType", pfType, &b_pfType);
   //fChain->SetBranchAddress("pfCandPt", pfCandPt, &b_pfCandPt);
   //fChain->SetBranchAddress("pfSumEcal", pfSumEcal, &b_pfSumEcal);
   //fChain->SetBranchAddress("pfSumHcal", pfSumHcal, &b_pfSumHcal);
   fChain->SetBranchAddress("trkStatus", trkStatus, &b_trkStatus);
   fChain->SetBranchAddress("nParticle", &nParticle, &b_nParticle);
   //fChain->SetBranchAddress("pStatus", pStatus, &b_pStatus);
   //fChain->SetBranchAddress("pPId", pPId, &b_pPId);
   fChain->SetBranchAddress("pEta", pEta, &b_pEta);
   fChain->SetBranchAddress("pPhi", pPhi, &b_pPhi);
   fChain->SetBranchAddress("pPt", pPt, &b_pPt);
   //fChain->SetBranchAddress("pAcc", pAcc, &b_pAcc);
   //fChain->SetBranchAddress("pAccPair", pAccPair, &b_pAccPair);
   //fChain->SetBranchAddress("pNRec", pNRec, &b_pNRec);
   //fChain->SetBranchAddress("pNHit", pNHit, &b_pNHit);
   fChain->SetBranchAddress("mtrkPt", mtrkPt, &b_mtrkPt);
   fChain->SetBranchAddress("mtrkPtError", mtrkPtError, &b_mtrkPtError);
   //fChain->SetBranchAddress("mtrkNHit", mtrkNHit, &b_mtrkNHit);
   //fChain->SetBranchAddress("mtrkNlayer", mtrkNlayer, &b_mtrkNlayer);
   //fChain->SetBranchAddress("mtrkNlayer3D", mtrkNlayer3D, &b_mtrkNlayer3D);
   fChain->SetBranchAddress("mtrkQual", mtrkQual, &b_mtrkQual);
   //fChain->SetBranchAddress("mtrkChi2", mtrkChi2, &b_mtrkChi2);
   //fChain->SetBranchAddress("mtrkNdof", mtrkNdof, &b_mtrkNdof);
   fChain->SetBranchAddress("mtrkDz1", mtrkDz1, &b_mtrkDz1);
   fChain->SetBranchAddress("mtrkDzError1", mtrkDzError1, &b_mtrkDzError1);
   fChain->SetBranchAddress("mtrkDxy1", mtrkDxy1, &b_mtrkDxy1);
   fChain->SetBranchAddress("mtrkDxyError1", mtrkDxyError1, &b_mtrkDxyError1);
   //fChain->SetBranchAddress("mtrkAlgo", mtrkAlgo, &b_mtrkAlgo);
   //fChain->SetBranchAddress("mtrkPfType", mtrkPfType, &b_mtrkPfType);
   //fChain->SetBranchAddress("mtrkPfCandPt", mtrkPfCandPt, &b_mtrkPfCandPt);
   //fChain->SetBranchAddress("mtrkPfSumEcal", mtrkPfSumEcal, &b_mtrkPfSumEcal);
   //fChain->SetBranchAddress("mtrkPfSumHcal", mtrkPfSumHcal, &b_mtrkPfSumHcal);
   Notify();
}

Bool_t trackTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trackTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trackTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void trackTree::Close()
{
f->Close();
}
#endif // #ifdef trackTree_cxx
