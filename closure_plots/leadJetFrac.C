#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TLatex.h"

void leadJetFrac()
{
  TH1::SetDefaultSumw2();

  const int nbins = 23;
  double bins[nbins] = {0.5,1,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};

  TFile * inf = new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/Correction_Vs3Calo_ntuple_dijet3.root","read");
  TNtuple * gen = (TNtuple*) inf->Get("nt_particle");
  TNtuple * reco =(TNtuple*) inf->Get("nt_track");

  TCanvas * c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,1,0,0);
  c1->cd(1);

  TH1D * gen_lead = new TH1D("gen_lead",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * gen_net = new TH1D("gen_net",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  
  gen->Draw("pt>>gen_lead","r_lead<0.3");
  gen->Draw("pt>>gen_net","r_lead<0.3 || r_sublead<0.3");
  gen_lead->Divide(gen_net);

  TH1D * reco_lead = new TH1D("reco_lead",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * reco_net = new TH1D("reco_net",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  
  reco->Draw("pt>>reco_lead","r_lead<0.3");
  reco->Draw("pt>>reco_net","r_lead<0.3 || r_sublead<0.3");
  reco_lead->Divide(reco_net);

  TH1D * reco_lead_corr = new TH1D("reco_lead_corr",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * reco_net_corr = new TH1D("reco_net_corr",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);

  reco->Draw("((1-fake)/eff)*pt>>reco_lead_corr","r_lead<0.3");
  reco->Draw("((1-fake)/eff)*pt>>reco_net_corr","r_lead<0.3 || r_sublead<0.3");
  reco_lead_corr->Divide(reco_net_corr);

  gen_lead->SetMarkerSize(0.8);
  gen_lead->SetLineWidth(1);
  gen_lead->SetMaximum(1);
  gen_lead->SetMinimum(0);
  gen_lead->GetXaxis()->SetTitle("p_{t}^{track}");
  gen_lead->GetYaxis()->SetTitle("(lead cone)/(dijet cone)");

  reco_lead->SetMarkerSize(0.8);
  reco_lead->SetLineWidth(1);
  reco_lead->SetMarkerColor(kRed+1);
  reco_lead->SetLineColor(kRed+1);

  reco_lead_corr->SetMarkerSize(0.8);
  reco_lead_corr->SetLineWidth(1);
  reco_lead_corr->SetMarkerColor(kBlue);
  reco_lead_corr->SetLineColor(kBlue);

  gen_lead->Draw("p e");
  reco_lead->Draw("same p e");
  reco_lead_corr->Draw("same p e");

  c1->cd(1);
  TLatex * lat = new TLatex(0.5,0.5," ");
  lat->DrawLatex(5,0.1,"Fragmentation JEC");
  TLegend * leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(gen_lead,"gen particles","p e");
  leg->AddEntry(reco_lead,"uncorr. tracks","p e");
  leg->AddEntry(reco_lead_corr,"corr. tracks","p e");
  leg->Draw("same");

  c1->cd(2);

  TFile * inf_2 = new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/Correction_Vs3Calo_ntuple_dijet_noFFJEC.root","read");
  TNtuple * gen_2 = (TNtuple*) inf_2->Get("nt_particle");
  TNtuple * reco_2 =(TNtuple*) inf_2->Get("nt_track");
 
  TH1D * gen_lead_2 = new TH1D("gen_lead_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * gen_net_2 = new TH1D("gen_net_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);

  gen_2->Draw("pt>>gen_lead_2","r_lead<0.3");
  gen_2->Draw("pt>>gen_net_2","r_lead<0.3 || r_sublead<0.3");
  gen_lead_2->Divide(gen_net_2);

  TH1D * reco_lead_2 = new TH1D("reco_lead_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * reco_net_2 = new TH1D("reco_net_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);

  reco_2->Draw("pt>>reco_lead_2","r_lead<0.3");
  reco_2->Draw("pt>>reco_net_2","r_lead<0.3 || r_sublead<0.3");
  reco_lead_2->Divide(reco_net_2);

  TH1D * reco_lead_corr_2 = new TH1D("reco_lead_corr_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);
  TH1D * reco_net_corr_2 = new TH1D("reco_net_corr_2",":p_{t}^{track}:(lead cone)/(dijet cone)",nbins-1,bins);

  reco_2->Draw("((1-fake)/eff)*pt>>reco_lead_corr_2","r_lead<0.3");
  reco_2->Draw("((1-fake)/eff)*pt>>reco_net_corr_2","r_lead<0.3 || r_sublead<0.3");
  reco_lead_corr_2->Divide(reco_net_corr_2);

  gen_lead_2->SetMarkerSize(0.8);
  gen_lead_2->SetLineWidth(1);
  gen_lead_2->SetMaximum(1);
  gen_lead_2->SetMinimum(0);
  gen_lead_2->GetXaxis()->SetTitle("p_{t}^{track}");
  gen_lead_2->GetYaxis()->SetTitle("(lead cone)/(dijet cone)");

  reco_lead_2->SetMarkerSize(0.8);
  reco_lead_2->SetLineWidth(1);
  reco_lead_2->SetMarkerColor(kRed+1);
  reco_lead_2->SetLineColor(kRed+1);

  reco_lead_corr_2->SetMarkerSize(0.8);
  reco_lead_corr_2->SetLineWidth(1);
  reco_lead_corr_2->SetMarkerColor(kBlue);
  reco_lead_corr_2->SetLineColor(kBlue);

  gen_lead_2->Draw("p e");
  reco_lead_2->Draw("same p e");
  reco_lead_corr_2->Draw("same p e");
  lat->DrawLatex(5,0.1,"No Fragmentation JEC");
}
