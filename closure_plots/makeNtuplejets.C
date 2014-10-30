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

//weightings for samples produced on 09_21_2014
float pthatWeight[7] = {0,0,0.000281494,5.95379e-05,5.93536e-05,5.81032e-05,6.11753e-05};
  float vertexShift = 0.436781;

 TString directory="/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/";
 

 const char* infname[7];
 infname[0] = "/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0";
 infname[1] = "/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0";
 infname[2] = "/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0";
 infname[3] = "/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0";
 infname[4] = "/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0";
 infname[5] = "/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0";
 infname[6] = "/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0";

 //full sample would be 350000,150000
 const int nevents[7] = {0,0,120000,0,0,0,0};
 
 trackTree * ftrk[7];
 HiTree * fhi[7];
 t * fjet[7];
 TFile * evtSelFile[7];
 TTree * evtSel[7];
 int pcoll[7];
 for(int ifile=2; ifile<3; ifile++){
   ftrk[ifile] = new trackTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile]));
   evtSelFile[ifile] = new TFile(Form("%s/%s.root",directory.Data(),infname[ifile]),"read");
   evtSel[ifile] = (TTree*) evtSelFile[ifile]->Get("skimanalysis/HltTree");
   evtSel[ifile]->SetBranchAddress("pcollisionEventSelection", &pcoll[ifile]);
  }

 TFile * centWeightsFile = new TFile("centrality_weights_MB.root","read");
 TH1F * centWeights = new TH1F("centWeight","centWeight",100,0,200);
 centWeights = (TH1F*)centWeightsFile->Get("centrality_weight");
 
 //pt bins for track efficiency correction 
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
   f_eff[ipt]= new TFile(Form("../final_hists_Vs3Calo/eff_pt%d_%d_cent%d_%d.root",(int)(100*ptmin_eff[ipt]),(int)(100*ptmax_eff[ipt]),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
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
   f_fake[ipt]= new TFile(Form("../final_hists_Vs3Calo/fake_pt%d_%d_cent%d_%d.root",(int)(100*ptmin_fake[ipt]),(int)(100*ptmax_fake[ipt]),(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt])));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }

 //output file and tree
 TFile *outf= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/Correction_Vs3Calo_ntuple.root","recreate");
 
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:cent:eff:cent_weight:pthat_weight:weight:pt1:pt2:dphi:asym:eta1:eta2:phi1:phi2:r_lead:r_sublead:isLeadClosest:isSubleadClosest";
 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:cent:eff:trkfake:fake:cent_weight:pthat_weight:weight:pt1:pt2:dphi:asym:eta1:eta2:phi1:phi2:r_lead:r_sublead:isLeadClosest:isSubleadClosest";
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());


 //loop over events

 for(int ifile=2; ifile<3; ifile++){
 std::cout<<ifile<<std::endl;
 int nentries = ftrk[ifile]->GetEntriesFast();
for(int jentry=0;jentry<nevents[ifile];jentry++){
 //for(int jentry=0;jentry<5000;jentry++){
  if((jentry%10000)==0) std::cout<<jentry<<"/"<<nevents[ifile]<<std::endl;
  ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;
 
  if(fabs(vz-vertexShift)>15 || !(pcoll[ifile])) continue;

  float weight = 0;
  float pthat_weight = 0;
  float cent_weight = 0;

  if(fjet[ifile]->pthat <50)       pthat_weight = pthatWeight[0];
  else if(fjet[ifile]->pthat <80)  pthat_weight = pthatWeight[1];
  else if(fjet[ifile]->pthat <120) pthat_weight = pthatWeight[2];
  else if(fjet[ifile]->pthat <220) pthat_weight = pthatWeight[3];
  else if(fjet[ifile]->pthat <280) pthat_weight = pthatWeight[4];
  else if(fjet[ifile]->pthat <370) pthat_weight = pthatWeight[5];
  else                             pthat_weight = pthatWeight[6];

//pthat weight for MBi
//delete for non-MB

  cent_weight = centWeights->GetBinContent(centWeights->FindBin(cent));
  weight = pthat_weight*cent_weight;

  float pt1=-99;
  float phi1=-99;
  float eta1=-99;
  float refpt1=-99;
  float refeta1=-99;
  float refphi1=-99;
  float matchedpt1=-99;
  float matchedR1=-99;
  float trackMax1=-99;
  float pt2=-99;
  float phi2=-99;
  float eta2=-99;
  float refpt2=-99;
  float refphi2=-99;
  float refeta2=-99;
  float matchedpt2=-99;
  float matchedR2=-99;
  float trackMax2=-99;
  float pt3=-99;
  float phi3=-99;
  float eta3=-99;
  float refpt3=-99;
  float refeta3=-99;
  float refphi3=-99;
  float matchedpt3=-99;
  float matchedR3=-99;
  float trackMax3=-99;
  float dphi=-99;
  float ptratio=-99;
  float asym = -1;

std::vector<std::pair<float, std::pair<float,std::pair<float, std::pair<float,std::pair<float,std::pair<float,std::pair<float,std::pair<float,float> > > > > > > > > jets;
  int njet=0;
  for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){

   if(fabs(fjet[ifile]->jteta[ijet])>2) continue;
   jets.push_back(std::make_pair(fjet[ifile]->jtpt[ijet],std::make_pair(fjet[ifile]->jteta[ijet], std::make_pair(fjet[ifile]->jtphi[ijet], std::make_pair(fjet[ifile]->refpt[ijet],std::make_pair(fjet[ifile]->refeta[ijet],std::make_pair(fjet[ifile]->refphi[ijet],std::make_pair(fjet[ifile]->matchedPt[ijet],std::make_pair(fjet[ifile]->matchedR[ijet],fjet[ifile]->trackMax[ijet])))))))));
   njet++;
  }

  std::sort(jets.begin(),jets.end());

//cut to study jet effects, removing no-jet events
//  if(njet == 0) continue;

  if(njet>0){
   pt1=       jets[njet-1].first;
   eta1=      jets[njet-1].second.first;
   phi1=      jets[njet-1].second.second.first;
   refpt1=    jets[njet-1].second.second.second.first;
   refeta1=   jets[njet-1].second.second.second.second.first;
   refphi1=   jets[njet-1].second.second.second.second.second.first;
   matchedpt1=jets[njet-1].second.second.second.second.second.second.first;
   matchedR1= jets[njet-1].second.second.second.second.second.second.second.first;
   trackMax1= jets[njet-1].second.second.second.second.second.second.second.second;
if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second.first;
    refpt2=jets[njet-2].second.second.second.first;
    refeta2=jets[njet-2].second.second.second.second.first;
    refphi2=jets[njet-2].second.second.second.second.second.first;
    matchedpt2=jets[njet-2].second.second.second.second.second.second.first;
    matchedR2=jets[njet-2].second.second.second.second.second.second.second.first;
    trackMax2=jets[njet-2].second.second.second.second.second.second.second.second;
    dphi=acos(cos(phi1-phi2));
    ptratio=pt2/pt1;
    asym = (pt1-pt2)/(pt1+pt2);
    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second.first;
     refpt3=jets[njet-3].second.second.second.first;
     refeta3=jets[njet-3].second.second.second.second.first;
     refphi3=jets[njet-3].second.second.second.second.second.first;
     matchedpt3=jets[njet-3].second.second.second.second.second.second.first;
     matchedR3=jets[njet-3].second.second.second.second.second.second.second.first;
     trackMax3=jets[njet-3].second.second.second.second.second.second.second.second;
    }
   }
  }


  //loop over tracks
  for(int itrk=0;itrk<ftrk[ifile]->nParticle;itrk++){

   float trackselect=(ftrk[ifile]->mtrkQual[itrk] && fabs(ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk])<3 && (ftrk[ifile]->mtrkPtError[itrk]/ftrk[ifile]->mtrkPt[itrk])<0.1);
   float eta=ftrk[ifile]->pEta[itrk];

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

    float isLeadClosest = 0;
    float isSubleadClosest = 0;
    float r_lead    = sqrt(pow(eta-eta1,2)+pow(acos(cos(phi-phi1)),2));
    if(r_lead == rmin) isLeadClosest = 1;
    float r_sublead = sqrt(pow(eta-eta2,2)+pow(acos(cos(phi-phi2)),2));
    if(r_sublead == rmin) isSubleadClosest = 1;   

  //cut for high R_lead or R_sublead so I can make a large ntuple
  //  if(isLeadClosest == 0 && isSubleadClosest == 0) continue;
  //  if(rmin < 1.6) continue;
    


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
  
   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff,cent_weight,pthat_weight,weight,pt1,pt2,dphi,asym,eta1,eta2,phi1,phi2,r_lead,r_sublead,isLeadClosest,isSubleadClosest};

   nt_particle->Fill(entry);

  }
  
  for(int itrk=0;itrk<ftrk[ifile]->nTrk;itrk++){

   float trackselect=(ftrk[ifile]->highPurity[itrk] && fabs(ftrk[ifile]->trkDxy1[itrk]/ftrk[ifile]->trkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->trkDz1[itrk]/ftrk[ifile]->trkDzError1[itrk])<3 && (ftrk[ifile]->trkPtError[itrk]/ftrk[ifile]->trkPt[itrk])<0.1);
   float eta=ftrk[ifile]->trkEta[itrk];

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


    float isLeadClosest = 0;
    float isSubleadClosest = 0;
    float r_lead    = sqrt(pow(eta-eta1,2)+pow(acos(cos(phi-phi1)),2));
    if(r_lead == rmin) isLeadClosest = 1;
    float r_sublead = sqrt(pow(eta-eta2,2)+pow(acos(cos(phi-phi2)),2));
    if(r_sublead == rmin) isSubleadClosest = 1;

  //cut for high R_lead or R_sublead so I can make a large ntuple
  //    if(isLeadClosest == 0 && isSubleadClosest == 0) continue;
  //    if(rmin < 1.6) continue;


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
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,cent,eff,trkfake,fake,cent_weight,pthat_weight,weight,pt1,pt2,dphi,asym,eta1,eta2,phi1,phi2,r_lead,r_sublead,isLeadClosest,isSubleadClosest};
   nt_track->Fill(entry);
  }
 }
}
 
  //nt_track->Write();
 // nt_particle->Write();
    outf->Write();
    outf->Close();
}
