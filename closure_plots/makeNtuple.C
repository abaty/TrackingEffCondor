#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "../ntupler/trackTree.C"

void makeNtuple(){
 TH1D::SetDefaultSumw2();

//all statistics 
 //float pthatWeight[5] = {4.29284e-01,2.99974e-02,3.38946e-4,1.06172e-4,2.79631e-5};
//100k events
// float pthatWeight[5] = {0.429284,0.0299974,0.000949812,0.000232709,7.61038e-05};
 //float vertexShift = 0.501501;
//100k events high pthat

 float pthatWeight[5] = {0,0,0.000654317,0.000156607,5.07966e-05};
 float vertexShift = 0.406408;

//current weights are for 150k (using 20k for approximation)
 const int nevents = 20000;

 TString directory="/mnt/hadoop/cms/store/user/velicanu/";
 const char* infname[5];
 infname[0] = "/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0/0";
 infname[1] = "/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0/0";
 infname[2] = "/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0/0";
 infname[3] = "/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0/0";
 infname[4] = "/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0/0";

 trackTree * ftrk[5];
 HiTree * fhi[5];
 t * fjet[5];
 TFile * evtSelFile[5];
 TTree * evtSel[5];
 int pcoll[5];
 for(int ifile=0; ifile<5; ifile++){
   ftrk[ifile] = new trackTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile]));
   evtSelFile[ifile] = new TFile(Form("%s/%s.root",directory.Data(),infname[ifile]),"read");
   evtSel[ifile] = (TTree*) evtSelFile[ifile]->Get("skimanalysis/HltTree");
   evtSel[ifile]->SetBranchAddress("pcollisionEventSelection", &pcoll[ifile]);
  }

 TFile * centWeightsFile = new TFile("centrality_weights.root","read");
 TH1F * centWeights = new TH1F("centWeight","centWeight",100,0,200);
 centWeights = (TH1F*)centWeightsFile->Get("centrality_weight");
 
 //pt bins for track efficiency correction Yen-Jie
 int npt_fake=29; 
 double ptmin_fake[]={0.5 ,0.5 ,0.5 ,0.5 ,0.5 ,0.55 ,0.55 ,0.55 ,0.55 ,0.55 ,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8,1 ,1 ,1 ,1 ,1  ,3 ,3 ,3  ,8};
 double ptmax_fake[]={0.55,0.55,0.55,0.55,0.55,0.65 ,0.65 ,0.65 ,0.65 ,0.65 ,0.8 ,0.8 ,0.8 ,0.8 ,0.8 ,1  ,1  ,1  ,1  ,1  ,3 ,3 ,3 ,3 ,3  ,8 ,8 ,8  ,300};
 
 int cent_min_fake[]={0   ,20  ,40  ,60  ,100  ,0    ,20   ,40   ,60   ,100   ,0   ,20  ,40  ,60  ,100  ,0  ,20 ,40 ,60 ,100 ,0 ,20,40,60,100 ,0 ,20,40 ,0};
 int cent_max_fake[]={20  ,40  ,60  ,100  ,200 ,20   ,40   ,60   ,100   ,200  ,20  ,40  ,60  ,100  ,200 ,20 ,40 ,60 ,100 ,200,20,40,60,100,200,20,40,200,200};
 
 int npt_eff=29; 
 double ptmin_eff[]={0.5 ,0.5 ,0.5 ,0.5 ,0.5 ,0.55 ,0.55 ,0.55 ,0.55 ,0.55 ,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8,1 ,1 ,1 ,1 ,1  ,3 ,3 ,3  ,8};
 double ptmax_eff[]={0.55,0.55,0.55,0.55,0.55,0.65 ,0.65 ,0.65 ,0.65 ,0.65 ,0.8 ,0.8 ,0.8 ,0.8 ,0.8 ,1  ,1  ,1  ,1  ,1  ,3 ,3 ,3 ,3 ,3  ,8 ,8 ,8  ,300};
 
 int cent_min[]={0   ,20  ,40  ,60  ,100  ,0    ,20   ,40   ,60   ,100   ,0   ,20  ,40  ,60  ,100  ,0  ,20 ,40 ,60 ,100 ,0 ,20,40,60,100 ,0 ,20,40 ,0};
 int cent_max[]={20  ,40  ,60  ,100  ,200 ,20   ,40   ,60   ,100   ,200  ,20  ,40  ,60  ,100  ,200 ,20 ,40 ,60 ,100 ,200,20,40,60,100,200,20,40,200,200};
  cout<<0<<endl;

 //getting histograms for track efficiency correction 
 TFile *f_eff[npt_eff];
 TProfile *p_eff_cent[npt_eff]; 
 TProfile2D *p_eff_accept[npt_eff]; 
 TProfile *p_eff_pt[npt_eff]; 
 TProfile *p_eff_rmin[npt_eff]; 
 for(int ipt=0; ipt<npt_eff;ipt++){
   f_eff[ipt]= new TFile(Form("../final_hists/eff_pt%d_%d_cent%d_%d.root",(int)(100*ptmin_eff[ipt]),(int)(100*ptmax_eff[ipt]),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt_fake];
 TProfile *p_fake_cent[npt_fake]; 
 TProfile2D *p_fake_accept[npt_fake]; 
 TProfile *p_fake_pt[npt_fake]; 
 TProfile *p_fake_rmin[npt_fake]; 
 for(int ipt=0; ipt<npt_fake;ipt++){
   f_fake[ipt]= new TFile(Form("../final_hists/fake_pt%d_%d_cent%d_%d.root",(int)(100*ptmin_fake[ipt]),(int)(100*ptmax_fake[ipt]),(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt])));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }

 //output file and tree
 TFile *outf= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_pthatCombo_20k_dzstudy.root","recreate");
 
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:cent:eff:cent_weight:pthat_weight:weight:dz_dzerr:dxy_dxyerr:Nhit:vz:trkAlgo";
 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:cent:eff:trkfake:fake:cent_weight:pthat_weight:weight:dz_dzerr:dxy_dxyerr:Nhit:vz:trkAlgo";
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());






 //loop over events

 for(int ifile=2; ifile<5; ifile++){
 std::cout<<ifile<<std::endl;
 int nentries = ftrk[ifile]->GetEntriesFast();
 if(nevents<nentries) nentries = nevents; 
for(int jentry=0;jentry<nentries;jentry++){
 //for(int jentry=0;jentry<5000;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;
 
  if(fabs(vz-vertexShift)>15 || !(pcoll[ifile])) continue;
  vz = vz-vertexShift;

  float weight = 0;
  float pthat_weight = 0;
  float cent_weight = 0;

  if(fjet[ifile]->pthat <50)       pthat_weight = pthatWeight[0];
  else if(fjet[ifile]->pthat <80)  pthat_weight = pthatWeight[1];
  else if(fjet[ifile]->pthat <100) pthat_weight = pthatWeight[2];
  else if(fjet[ifile]->pthat <120) pthat_weight = pthatWeight[3];
  else                             pthat_weight = pthatWeight[4];

  cent_weight = centWeights->GetBinContent(centWeights->FindBin(cent));
  weight = pthat_weight*cent_weight;

  //loop over tracks
  for(int itrk=0;itrk<ftrk[ifile]->nParticle;itrk++){

   float trackselect=(ftrk[ifile]->mtrkQual[itrk] && fabs(ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk])<3 && (ftrk[ifile]->mtrkPtError[itrk]/ftrk[ifile]->mtrkPt[itrk])<0.1);
   float dz_cut=fabs(ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk]);
   float dxy_cut=fabs(ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk]);
   float Nhit = ftrk[ifile]->trkNHit[itrk];
   float eta=ftrk[ifile]->pEta[itrk];
   float trkAlgo=ftrk[ifile]->trkAlgo[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=ftrk[ifile]->pPt[itrk];
   float mpt=ftrk[ifile]->mtrkPt[itrk];
   float phi=ftrk[ifile]->pPhi[itrk];
   float rmin=199;
 
   for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
     if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
     float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
     if(r_reco<rmin)rmin=r_reco;
    }

   
   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 

   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   
   //fill in the output tree
  
   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff,cent_weight,pthat_weight,weight,dz_cut,dxy_cut,Nhit,vz,trkAlgo};

   nt_particle->Fill(entry);

  }
  
  for(int itrk=0;itrk<ftrk[ifile]->nTrk;itrk++){

   float trackselect=(ftrk[ifile]->highPurity[itrk] && fabs(ftrk[ifile]->trkDxy1[itrk]/ftrk[ifile]->trkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->trkDz1[itrk]/ftrk[ifile]->trkDzError1[itrk])<3 && (ftrk[ifile]->trkPtError[itrk]/ftrk[ifile]->trkPt[itrk])<0.1);
   float dz_cut=fabs(ftrk[ifile]->trkDz1[itrk]/ftrk[ifile]->trkDzError1[itrk]);
   float dxy_cut=fabs(ftrk[ifile]->trkDxy1[itrk]/ftrk[ifile]->trkDxyError1[itrk]);
   float Nhit = ftrk[ifile]->trkNHit[itrk];
   float eta=ftrk[ifile]->trkEta[itrk];
   float trkAlgo=ftrk[ifile]->trkAlgo[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker   
   
   float pt=ftrk[ifile]->trkPt[itrk];
   float phi=ftrk[ifile]->trkPhi[itrk];
   float trkfake=ftrk[ifile]->trkFake[itrk];
   float trackstatus=ftrk[ifile]->trkStatus[itrk];
   float rmin=199;

   //find rmin; 
     for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
     if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
     float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
     if(r_reco<rmin)rmin=r_reco;
    }

   //get efficiency and fake rate correction for the track Yen-Jie
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;
   
   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin)); 
     }     
   } 
   
   for(int ipt=0;ipt<npt_fake;ipt++){
    if(pt>=ptmin_fake[ipt] && pt<ptmax_fake[ipt] && cent>=cent_min_fake[ipt] && cent<cent_max_fake[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<=100) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }

   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   if(eff==0){
          //cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }

   //if(fake<0) fake=0;

   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,cent,eff,trkfake,fake,cent_weight,pthat_weight,weight,dz_cut,dxy_cut,Nhit,vz,trkAlgo};
   nt_track->Fill(entry);
  }
 }
}
 
  //nt_track->Write();
 // nt_particle->Write();
    outf->Write();
    outf->Close();
}
