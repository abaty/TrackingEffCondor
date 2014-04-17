#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"

void plot_efficiency_cent(int nstep_cent=2,int nstep_accept=2,int nstep_pt=2,int nstep_rmin=2,double  bin_pt_min=0.4,double bin_pt_max=1,double bin_cent_min=0, double bin_cent_max=10,int nevents=45000,bool is_final=false){
TH1D::SetDefaultSumw2();
TH2D::SetDefaultSumw2(true);

// TFile * f = new TFile(Form("../workDir/tracking_ntuples/new/track_ntuple_Dijet100_HydjetDrum_v27_mergedV1_%devts_cent_%d_accept_%d_pt_%d_ptmin%d_ptmax%d.root",nevents,nstep_cent,nstep_accept,nstep_pt,bin_pt_min,bin_pt_max));
TFile * f = new TFile(Form("track_ntuple_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_centmin%d_centmax%d.root",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int) bin_cent_max));
TTree * nt = (TTree*)f->Get("nt_track");

//##########################################Initial Efficiencies without Correction###########################################################################
///////////acceptance dependent efficiency/////////
TProfile2D * p_eta_phi = new TProfile2D(Form("p_eta_phi"),";#phi;#eta;",50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
nt->Draw(Form("(mpt>0 && trackselect):eta:phi>>p_eta_phi"),Form("weight*(abs(eta)<2.4 && pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

//////////pt dependent efficiency////////////

const int ny=20;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(20*log(10));
int maxbin=ny;
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix);
 if(x[ix]>100){
	x[ix]=bin_pt_max;
	maxbin=ix;
	break;
 }
 inix+=delta;
} 

TProfile * p_pt= new TProfile("p_pt",";p_{T}(GeV/c);efficiency",maxbin,x);
TH1D * pt_gen = new TH1D("pt_gen",";p_{T}(GeV/c);N_{part}",maxbin,x);
TH1D * pt_reco = new TH1D("pt_reco",";p_{T}(GeV/c);N_{part}",maxbin,x);
nt->Draw("pt>>pt_gen","weight");
nt->Draw("pt>>pt_reco","weight*(mpt>0)");
nt->Draw("(mpt>0  && trackselect):pt>>p_pt","weight");

pt_gen->SetMarkerStyle(25);
pt_gen->SetMinimum(pt_reco->GetMinimum());
pt_gen->SetMaximum(1.2*pt_gen->GetMaximum());

////////centrality dependent efficiency//////////

TProfile * p_cent= new TProfile("p_cent",";centrality bin;efficiency",100,0,200);
TH1D * cent_gen = new TH1D("cent_gen",";centrality bin;N_{part}",100,0,200);
TH1D * cent_reco = new TH1D("cent_reco",";centrality bin;N_{part}",100,0,200);
nt->Draw("cent>>cent_gen","weight");
nt->Draw("cent>>cent_reco","weight*(mpt>0)");
nt->Draw("(mpt>0 && trackselect):cent>>p_cent",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TH1D * eff_cent=(TH1D*)cent_reco->Clone("eff_cent");
eff_cent->Divide(cent_gen);

cent_gen->SetMarkerStyle(25);
cent_gen->SetMinimum(cent_gen->GetMinimum());
cent_gen->SetMaximum(1.2*cent_gen->GetMaximum());

//################################################rmin efficiency (calculated elsewhere, for bookkeeping here#################################################

const int n_rmin_bins=36;
double rmin_bins[n_rmin_bins+1] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3,100};

TProfile * p_rmin = new TProfile("p_rmin",";#phi;#eta;efficiency",n_rmin_bins,rmin_bins);
nt->Draw("eff_rmin:rmin_reco>>p_rmin",Form("weight*( abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

//###############################################Efficiencies after Correction###########################################################################
///////////////cent dependent///////////////////
TProfile * p_cent_corr= new TProfile("p_cent_corr",";centrality bin;efficiency",100,0,200);
nt->Draw("(1/eff)*(mpt>0 && trackselect):cent>>p_cent_corr",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

p_cent_corr->SetMaximum(1.1);
p_cent_corr->SetMinimum(0.9);

////////////pt dependent////////////////////////
TProfile * p_pt_corr= new TProfile("p_pt_corr",";p_{T}(GeV/c);efficiency",maxbin,x);
nt->Draw("(1/eff)*(mpt>0 && trackselect):pt>>p_pt_corr",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

p_pt_corr->SetMaximum(1.1);
p_pt_corr->SetMinimum(0.9);

//////////acceptance dependent///////////////////

TProfile2D * p_eta_phi_corr = new TProfile2D("p_eta_phi_corr",";#phi;#eta;",50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
nt->Draw("(1/eff)*(mpt>0 && trackselect):eta:phi>>p_eta_phi_corr",Form("weight*( abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

////////////rmin dependent/////////////////////////
TProfile * p_rmin_corr = new TProfile("p_rmin_corr",";r_{min};efficiency",n_rmin_bins,rmin_bins);
nt->Draw("(1/eff)*(mpt>0 && trackselect):rmin_reco>>p_rmin_corr",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

p_rmin_corr->SetMaximum(1.1);
p_rmin_corr->SetMinimum(0.9);

/////output file to be used in the ntupler in the next step/////
TFile *outf = new TFile(Form("eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max),"recreate");
p_cent->Write();
p_pt->Write();
p_eta_phi->Write();
p_cent_corr->Write();
p_pt_corr->Write();
p_eta_phi_corr->Write();
p_rmin_corr->Write();
p_rmin->Write();

outf->Close();

////overall efficiency histograms///////////////////////////////
TProfile * p_eff_cent = new TProfile("p_eff_cent",";centrality bin;efficiency",100,0,200);
nt->Draw("eff_cent:cent>>p_eff_cent",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile * p_eff_pt = new TProfile("p_eff_pt",";p_{T} bin;efficiency",maxbin,x);
nt->Draw("eff_pt:pt>>p_eff_pt",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile2D * p_eff_acceptance = new TProfile2D("p_eff_acceptance",";#phi;#eta;efficiency",50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
nt->Draw("eff_accept:eta:phi>>p_eff_acceptance",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile * p_eff_rmin = new TProfile("p_eff_rmin",";#phi;#eta;efficiency",n_rmin_bins,rmin_bins);
nt->Draw("eff_rmin:rmin_reco>>p_eff_rmin",Form("weight*(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");


TFile *f_efficiency;
 if(is_final){
 f_efficiency= new TFile(Form("eff_pt%d_%d_cent%d_%d_step_cent%daccept%dpt%drmin%d.root",(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max,nstep_cent,nstep_accept,nstep_pt,nstep_rmin),"recreate");

 p_eff_cent->SetMarkerStyle(20);
 p_eff_pt->SetMarkerStyle(20);
 p_eff_rmin->SetMarkerStyle(20);

 p_eff_cent->Write();
 p_eff_pt->Write();
 p_eff_acceptance->Write();
 p_eff_rmin->Write();

 f_efficiency->Close();
 f->Close();

//memory cleanup
 delete p_eff_cent;
 delete p_eff_pt;
 delete p_eff_acceptance;
 delete p_eff_rmin;
}
}
