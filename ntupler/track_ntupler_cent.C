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

void track_ntupler_cent(int nstep_cent=2,int nstep_accept=1,int nstep_pt=1,int nstep_rmin=1,double bin_pt_min=8,double bin_pt_max=100,double bin_cent_min=0,double bin_cent_max=10,int *nevents=0){
 TH1D::SetDefaultSumw2();

 //converted to nb
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

 trackTree * ftrk[7];
 HiTree * fhi[7];
 t * fjet[7];
 TFile * evtSelFile[7];
 TTree * evtSel[7];
 int pcoll[7];
 for(int ifile=0; ifile<7; ifile++){
   ftrk[ifile] = new trackTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile]));
   evtSelFile[ifile] = new TFile(Form("%s/%s.root",directory.Data(),infname[ifile]),"read");
   evtSel[ifile] = (TTree*) evtSelFile[ifile]->Get("skimanalysis/HltTree");
   evtSel[ifile]->SetBranchAddress("pcollisionEventSelection", &pcoll[ifile]); 
  }

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

 TFile * centWeightsFile = new TFile("centrality_weights.root","read");
 TH1F * centWeights = new TH1F("centWeight","centWeight",100,0,200);
 centWeights = (TH1F*)centWeightsFile->Get("centrality_weight");

 TFile *outf= new TFile(Form("track_ntuple_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_centmin%d_centmax%d.root",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max),"recreate");

 std::string trackVars="pthat:pt:mpt:eta:phi:trackselect:cent:rmin_reco:eff_cent:eff_accept:eff_pt:eff_rmin:eff:weight:pthat_weight:cent_weight";

 
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());


//file loop here
//note we only use files 2-4 for now for stats reasons
  for(int ifile=2; ifile<4; ifile++){
  std::cout<<ifile<<std::endl;

//event loop 
  for(int jentry=0;jentry<nevents[ifile];jentry++){
  if((jentry%10000)==0) std::cout<<jentry<<"/"<<nevents[ifile]<< "   File:" << ifile <<std::endl;

  ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;

//vertexShift can be used to move the center of the vertex distribution if needed 
  if(cent*0.5<bin_cent_min || cent*0.5>=bin_cent_max || fabs(vz-vertexShift)>15 || !(pcoll[ifile])) continue;
 
  float eff_cent=1;

  for(int icent=0;icent<nstep_cent;icent++){
    eff_cent=eff_cent*p_eff_cent[icent]->GetBinContent(p_eff_cent[icent]->FindBin(cent));
  }

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

  cent_weight = centWeights->GetBinContent(centWeights->FindBin(cent));  
  weight = pthat_weight*cent_weight;

//particle loop here
  for(int itrk=0;itrk<ftrk[ifile]->nParticle;itrk++){
   
   float pt=ftrk[ifile]->pPt[itrk];
   if(pt<bin_pt_min || pt>bin_pt_max) continue;

   float eta=ftrk[ifile]->pEta[itrk];
   if(fabs(eta)>2.4) continue;

   float trackselect=(ftrk[ifile]->mtrkQual[itrk] && fabs((ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk]))<3.0 && fabs((ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk]))<3 && (ftrk[ifile]->mtrkPtError[itrk]/ftrk[ifile]->mtrkPt[itrk])<0.1);

   float mpt=ftrk[ifile]->mtrkPt[itrk];
   float phi=ftrk[ifile]->pPhi[itrk];
   
   float eff_accept=1;
   float eff_pt=1;
   float eff_rmin=1;
   float rmin_reco=199;
  
   for(int iaccept=0;iaccept<nstep_accept;iaccept++){
    eff_accept=eff_accept*p_eff_accept[iaccept]->GetBinContent(p_eff_accept[iaccept]->GetXaxis()->FindBin(phi),p_eff_accept[iaccept]->GetYaxis()->FindBin(eta)); 
   }
  
   for(int ipt=0;ipt<nstep_pt;ipt++){
    eff_pt=eff_pt*p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->GetXaxis()->FindBin(pt));    
   } 
   
   for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
    if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
    float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
    if(r_reco<rmin_reco)rmin_reco=r_reco;
   }
   
   for(int irmin=0;irmin<nstep_rmin;irmin++){
     if(rmin_reco<=100)eff_rmin=eff_rmin*p_eff_rmin[irmin]->GetBinContent(p_eff_rmin[irmin]->GetXaxis()->FindBin(rmin_reco));
   }
   
   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;

   float entry[]={fjet[ifile]->pthat,pt,mpt,eta,phi,trackselect,cent,rmin_reco,eff_cent,eff_accept,eff_pt,eff_rmin,eff,weight,pthat_weight,cent_weight};
   nt_track->Fill(entry);
  }
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
  for(int ifile=0; ifile<7; ifile++){		
    fjet[ifile]->Close();
    fhi[ifile]->Close();
    ftrk[ifile]->Close();
    evtSelFile[ifile]->Close();
  }
  centWeightsFile->Close();
}
