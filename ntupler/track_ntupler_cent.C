#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
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
#include "trackTree.C"

void track_ntupler_cent(int nstep_cent=2,int nstep_accept=1,int nstep_pt=1,int nstep_rmin=1,double bin_pt_min=8,double bin_pt_max=100,double bin_cent_min=0,double bin_cent_max=10,int nevents=593463){
 TH1D::SetDefaultSumw2();
 double R=0.3;
 
 //TString directory="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/";
 //TString infname="Dijet100_HydjetDrum_v27_mergedV1";
 TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet50/Track8_Jet19_STARTHI53_LV1/merged/";
 TString infname="HiForest_Pythia_Hydjet_Jet50_Track8_Jet19_STARTHI53_LV1_merged_forest_0"; 

 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()));
 
 TFile *f_eff_cent[nstep_cent];
 TProfile *p_eff_cent[nstep_cent]; 
 for(int icent=0; icent<nstep_cent;icent++){
   f_eff_cent[icent]= new TFile(Form("eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",icent, icent, icent,icent,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max));
   p_eff_cent[icent]=(TProfile*)f_eff_cent[icent]->Get("p_cent_corr");
 }

 TFile *f_eff_accept[nstep_accept];
 TProfile2D * p_eff_accept[nstep_accept];
 for(int iaccept=0;iaccept<nstep_accept;iaccept++){
   f_eff_accept[iaccept]= new TFile(Form("eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",iaccept+1, iaccept, iaccept,iaccept,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max));
   p_eff_accept[iaccept] = (TProfile2D*)f_eff_accept[iaccept]->Get("p_eta_phi_corr");
 }
 

 TFile *f_eff_pt[nstep_pt];
 TProfile * p_eff_pt[nstep_pt];
 
 for(int ipt=0;ipt<nstep_pt;ipt++){
   f_eff_pt[ipt]= new TFile(Form("eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",ipt+1, ipt+1, ipt,ipt,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max));
   p_eff_pt[ipt] = (TProfile*)f_eff_pt[ipt]->Get("p_pt_corr");
 }
 
  
 TFile *f_eff_rmin[nstep_rmin];
 TProfile * p_eff_rmin[nstep_rmin];
 
 for(int irmin=0;irmin<nstep_rmin;irmin++){
   f_eff_rmin[irmin]= new TFile(Form("eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",irmin+1, irmin+1, irmin+1,irmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max));
   p_eff_rmin[irmin] = (TProfile*)f_eff_rmin[irmin]->Get("p_rmin_corr");
 }
 
 
 int nentries = ftrk->GetEntriesFast();

 TFile *outf= new TFile(Form("track_ntuple_%s_%devts_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_centmin%d_centmax%d.root",infname.Data(),nevents,nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max),"recreate");
 // TFile *outf= new TFile(Form("../ntuples/track_ntuple_%s_%devts_cent_%d_accept_%d_pt_%d_ptmin%d_ptmax%d.root",infname.data(),nevents,nstep_cent,nstep_accept,nstep_pt,(int)bin_pt_min,(int)bin_pt_max),"recreate");

 // string trackVars="pt:mpt:eta:phi:trackselect:cent:vx:vy:vz:pt1:pt2:dphi:ptratio:incone1:incone2";
 std::string trackVars="pthat:pt:mpt:eta:phi:trackselect:cent:pt1:pt2:pt3:eta1:eta2:eta3:phi1:phi2:phi3:dphi:incone1:incone2:incone3:refpt1:refpt2:refpt3:refeta1:refeta2:refeta3:refphi1:refphi2:refphi3:matchedpt1:matchedpt2:matchedpt3:matchedR1:matchedR2:matchedR3:trackMax1:trackMax2:trackMax3:incone:rmin_reco:rmin_gen:eff_cent:eff_accept:eff_pt:eff_rmin:eff";
 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 // for(int jentry=0;jentry<nentries;jentry++){
  for(int jentry=0;jentry<nevents;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);

  float cent=fhi->hiBin;
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
  //float ptratio=-99;
  float eff_cent=1;
  refpt3=-99;

  if(cent*0.5<bin_cent_min || cent*0.5>=bin_cent_max) continue;
  for(int icent=0;icent<nstep_cent;icent++){
    eff_cent=eff_cent*p_eff_cent[icent]->GetBinContent(p_eff_cent[icent]->FindBin(cent));
  }
 
  std::vector<std::pair<double, std::pair<double,std::pair<double, std::pair<double,std::pair<double,std::pair<double,std::pair<double,std::pair<double,double> > > > > > > > > jets;
  int njet=0;
  for(int ijet=0;ijet<fjet->nref;ijet++){

   if(fabs(fjet->jteta[ijet])>2) continue;
   jets.push_back(std::make_pair(fjet->jtpt[ijet],std::make_pair(fjet->jteta[ijet], std::make_pair(fjet->jtphi[ijet], std::make_pair(fjet->refpt[ijet],std::make_pair(fjet->refeta[ijet],std::make_pair(fjet->refphi[ijet],std::make_pair(fjet->matchedPt[ijet],std::make_pair(fjet->matchedR[ijet],fjet->trackMax[ijet])))))))));
   njet++;

  }

  std::sort(jets.begin(),jets.end());
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
    //ptratio=pt2/pt1;
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

  for(int itrk=0;itrk<ftrk->nParticle;itrk++){
   float trackselect=(ftrk->mtrkQual[itrk] && fabs((ftrk->mtrkDxy1[itrk]/ftrk->mtrkDxyError1[itrk]))<3.0 && fabs((ftrk->mtrkDz1[itrk]/ftrk->mtrkDzError1[itrk]))<3 && (ftrk->mtrkPtError[itrk]/ftrk->mtrkPt[itrk])<0.1);
   
   float pt=ftrk->pPt[itrk];
   if(pt<bin_pt_min || pt>bin_pt_max) continue;
   float mpt=ftrk->mtrkPt[itrk];
   float eta=ftrk->pEta[itrk];
   if(fabs(eta)>2.4) continue;

   float phi=ftrk->pPhi[itrk];
   float incone1=0;
   float incone2=0;
   float incone3=0;
   float incone=0;
   
   float eff_accept=1;
   float eff_pt=1;
   float eff_rmin=1;
   float rmin_reco=99;
   float rmin_gen=99;

  
   for(int iaccept=0;iaccept<nstep_accept;iaccept++){
    eff_accept=eff_accept*p_eff_accept[iaccept]->GetBinContent(p_eff_accept[iaccept]->GetXaxis()->FindBin(phi),p_eff_accept[iaccept]->GetYaxis()->FindBin(eta)); 
   }
  
   for(int ipt=0;ipt<nstep_pt;ipt++){
    eff_pt=eff_pt*p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->GetXaxis()->FindBin(pt));    
   } 
  
   if(sqrt(pow(eta-eta1,2)+pow(acos(cos(phi-phi1)),2))<R) incone1=1;
   if(sqrt(pow(eta-eta2,2)+pow(acos(cos(phi-phi2)),2))<R) incone2=1;
   if(sqrt(pow(eta-eta3,2)+pow(acos(cos(phi-phi3)),2))<R) incone3=1;
   
   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    float r_gen=sqrt(pow(eta-fjet->refeta[ijet],2)+pow(acos(cos(phi-fjet->refphi[ijet])),2));
    if(r_reco<R) incone=1;
    if(r_reco<rmin_reco)rmin_reco=r_reco;
    if(r_gen<rmin_gen)rmin_gen=r_gen;
   }
   
   for(int irmin=0;irmin<nstep_rmin;irmin++){
    if(rmin_reco<5)eff_rmin=eff_rmin*p_eff_rmin[irmin]->GetBinContent(p_eff_rmin[irmin]->GetXaxis()->FindBin(rmin_reco));
   }
   
   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;

   float entry[]={fjet->pthat,pt,mpt,eta,phi,trackselect,cent,pt1,pt2,pt3,eta1,eta2,eta3,phi1,phi2,phi3,dphi,incone1,incone2,incone3,refpt1,refpt2,refpt3,refeta1,refeta2,refeta3,refphi1,refphi2,refphi3,matchedpt1,matchedpt2,matchedpt3,matchedR1,matchedR2,matchedR3,trackMax1,trackMax2,trackMax3,incone,rmin_reco,rmin_gen,eff_cent,eff_accept,eff_pt,eff_rmin,eff};
   nt_track->Fill(entry);
  }
 }
 
  nt_track->Write();
  outf->Close();

//some memory cleanup
  for(int icent=0; icent<nstep_cent;icent++){
  	f_eff_cent[icent]->Close();	
	}
  for(int iaccept=0; iaccept<nstep_accept;iaccept++){
        f_eff_accept[iaccept]->Close();
	}	
  for(int ipt=0; ipt<nstep_pt;ipt++){
        f_eff_pt[ipt]->Close();
	}	
  for(int irmin=0; irmin<nstep_rmin;irmin++){
        f_eff_rmin[irmin]->Close();
	}		
  fjet->Close();
  fhi->Close();
  ftrk->Close();
}
