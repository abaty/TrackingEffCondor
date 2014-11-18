#include <iostream>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

struct fragmentation_JEC
{
 private:
  static const int ncent=4;
  int radius;
  int ntrkmax;
  int cent_max[ncent];
  int cent_min[ncent];
  double PF_pt_cut;
  double PF_eta_cut;
  bool do_PbPb;
  bool do_pp_tracking;
  TString algo_corr; 
  TH2D* correction_matrix[ncent];
  TFile* correction_file;
 public:
  void reset()
  { 
   for(int icent=0;icent<ncent;icent++){
    correction_matrix[icent]=NULL;
   }
   correction_file=NULL;
   ntrkmax=21;
   PF_eta_cut=2.4;
   cent_min[0]=0;
   cent_max[0]=cent_min[1]=20;
   cent_max[1]=cent_min[2]=60;
   cent_max[2]=cent_min[3]=100;
   cent_max[3]=200;
  }
  
  fragmentation_JEC(int radius=3, bool do_PbPb=1, bool do_pp_tracking=0, double PF_pt_cut=2)
  {
   reset();
   if(do_PbPb==1){
    do_pp_tracking=0;
   }
   this->do_PbPb=do_PbPb;
   this->radius=radius;
   this->PF_pt_cut=PF_pt_cut;
   this->do_pp_tracking=do_pp_tracking;
   if(PF_pt_cut==3) ntrkmax=26;
   else if(PF_pt_cut==2) ntrkmax=31;
   else if(PF_pt_cut==1) ntrkmax=41;
  }
  
  double get_corrected_pt(double jetpt, int ntrk, int cent=0)
  {
   double correction=1;
   
   int cent_bin=0;
   if(do_PbPb){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   if(ntrk<ntrkmax) correction=correction_matrix[cent_bin]->GetBinContent(correction_matrix[cent_bin]->GetXaxis()->FindBin(ntrk),correction_matrix[cent_bin]->GetYaxis()->FindBin(jetpt));
   else correction=correction_matrix[cent_bin]->GetBinContent(correction_matrix[cent_bin]->GetXaxis()->FindBin(ntrkmax),correction_matrix[cent_bin]->GetYaxis()->FindBin(jetpt));
   
   return (1/correction)*jetpt;
  }
  
  bool passes_PF_selection(double PF_pt, double PF_eta, double PF_phi, int PF_id, double jet_eta, double jet_phi)
  { 
   if(PF_pt<PF_pt_cut || fabs(PF_eta)>PF_eta_cut || PF_id!=1) return false;
   
   double r=sqrt(pow(jet_eta-PF_eta,2)+pow(acos(cos(jet_phi-PF_phi)),2));
   
   if(r<((double)(radius)*0.1)) return true;
   else return false;
  }
  
  void set_correction()
  {
   std::cout<<"setting correction"<<std::endl;
   if(do_PbPb){
    algo_corr=Form("akVs%dCalo",radius);
    correction_file = new TFile(Form("FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    for(int icent=0;icent<ncent;icent++){
	 correction_matrix[icent]=(TH2D*)correction_file->Get(Form("pNtrk_pt%d",icent));
    } 
   }else{
  	algo_corr=Form("ak%dCalo",radius);
    if(do_pp_tracking){
     correction_file = new TFile(Form("FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     correction_matrix[0]=(TH2D*)correction_file->Get("pNtrk_pt");    
    }else{
     correction_file = new TFile(Form("FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     correction_matrix[0]=(TH2D*)correction_file->Get("pNtrk_pt");
    }
   }
  }
};
