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


void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth = 
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}
 
void fixAxes(TH1D* histo){
  histo->GetYaxis()->SetLabelSize(12);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleSize(20);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetTitleOffset(1.5);
  histo->GetYaxis()->CenterTitle();
  histo->GetXaxis()->SetLabelSize(12);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetXaxis()->SetTitleSize(20);
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleOffset(1.5);
}

void legFormat(TLegend* leg){
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0,"CMS Preliminary","");
  leg->AddEntry((TObject*)0,"PYTHIA+HYDJET","");
  leg->AddEntry((TObject*)0,"VsCalo jets, HI tracking","");
}

void drawClosure(TLine * l, TH1D * histo, int i=0){
histo->SetMaximum(1.1);
histo->SetMinimum(0.9);
if(i==0) histo->GetYaxis()->SetTitle("Efficiency Closure");
if(i==1)  histo->GetYaxis()->SetTitle("Fake Closure");
if(i==2)   histo->GetYaxis()->SetTitle("Total Closure");
histo->Draw();
fixAxes(histo);
l->Draw("same");
}

void pt_closure_test(int m = 0){

if(m<0 || m>4){
  m=0;
  std::cout << "invalid mode! Assuming pt mode." << std::endl;
}

const char* var[5] = {"pt","cent","rmin","eta","phi"};
const char* label[5] = {"p_{T}","cent","r_{min}","#eta","#phi"}; 

TH1D::SetDefaultSumw2();
TFile * f= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_pthatCombo_150kfull.root","read");
TTree * nt_track = (TTree*)f->Get("nt_track");
TTree * nt_particle = (TTree*)f->Get("nt_particle");

TLegend *leg = new TLegend(0.6,0.75,0.95,0.95);
leg->SetBorderSize(0);
leg->SetFillStyle(0);

double bin_pt_min=0.5;
double bin_pt_max=300;

const int ny=50;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(ny*log(10));
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix); 
 inix+=delta;
}
//Normalization factor

TH1D * h_weight = new TH1D("h_weight","h_weight",100,0,100);
nt_particle->Draw("weight>>h_weight","");
double weight_integral = h_weight->Integral(1,100);
std::cout << "normalization " <<weight_integral << std::endl;



//eff correction
TH1D * h_gen = new  TH1D("h_gen",";p_{T};Arbitrary Units",ny,x);
TH1D * h_gen_select = new  TH1D("h_gen_select",";p_{T};N_{evt}",ny,x);
TH1D * h_gen_matched_select_corr = new TH1D("h_gen_matched_select_corr",";p_{T};N_{evt}",ny,x);

nt_particle->Draw("pt>>h_gen","weight*(pt>0.5)");
nt_particle->Draw("pt>>h_gen_select","weight*(trackselect && pt>0.5)");
nt_particle->Draw("pt>>h_gen_matched_select_corr","(1/eff)*weight*(trackselect && pt>0.5)");

h_gen->SetMarkerColor(1);
h_gen->SetMarkerStyle(25);
h_gen->SetLineWidth(1);
h_gen_select->SetMarkerColor(1);
h_gen_matched_select_corr->SetMarkerColor(kRed);
h_gen_matched_select_corr->SetLineColor(kRed);

TLegend *leg2 = new TLegend(0.15,0.1,0.8,0.4);
legFormat(leg2);
leg2->AddEntry(h_gen,"gen","p");
leg2->AddEntry(h_gen_select,"matched gen","p");
leg2->AddEntry(h_gen_matched_select_corr,"corrected matched gen","p");

TCanvas * c2 = new TCanvas("c2","",600,600);
makeMultiPanelCanvas(c2,1,2,0.0,0.0,0.15,0.15,0.02);
c2->cd(1);
c2->cd(1)->SetLogx();
c2->cd(1)->SetLogy();

h_gen->Draw();
fixAxes(h_gen);
h_gen_select->Draw("same");
h_gen_matched_select_corr->Draw("same");
leg2->Draw("same");

TH1D * hgen_corr_rat = (TH1D*)h_gen_matched_select_corr->Clone("hgen_corr_rat");
hgen_corr_rat->Divide(h_gen);
c2->cd(2);
c2->cd(2)->SetLogx();

TLine * l = new TLine(0.5,1,300,1);
drawClosure(l,hgen_corr_rat,0);

c2->SaveAs("compare_gen_select_corr.png");
c2->SaveAs("compare_gen_select_corr.pdf");

//****************************************************************************
//fake correction
TH1D * h_reco = new  TH1D("h_reco",";p_{T};Arbitrary Units",ny,x);
TH1D * h_reco_matched = new  TH1D("h_reco_matched",";p_{T};N_{evt}",ny,x);
TH1D * h_reco_fakecorr = new TH1D("h_reco_fakecorr",";p_{T};N_{evt}",ny,x);
nt_track->Draw("pt>>h_reco","weight*(trackselect && pt>0.5)"); 
nt_track->Draw("pt>>h_reco_matched","weight*(trackselect && !trkfake && pt>0.5)"); 
nt_track->Draw("pt>>h_reco_fakecorr","(1-fake)*weight*(trackselect && pt>0.5)"); 

h_reco_matched->SetMarkerColor(1);
h_reco_matched->SetMarkerStyle(25);
h_reco_fakecorr->SetMarkerColor(kRed);
h_reco_fakecorr->SetLineColor(kRed);

TLegend *leg3 = new TLegend(0.15,0.1,0.8,0.4);
legFormat(leg3);
leg3->AddEntry(h_reco,"reco","p");
leg3->AddEntry(h_reco_matched,"matched reco","p");
leg3->AddEntry(h_reco_fakecorr,"fake corrected reco","p");

TCanvas * c3 = new TCanvas("c3","",600,600);
makeMultiPanelCanvas(c3,1,2,0.0,0.0,0.15,0.15,0.02);
c3->cd(1);
c3->cd(1)->SetLogx();
c3->cd(1)->SetLogy();
h_reco->Draw();
fixAxes(h_reco);

h_reco_matched->Draw("same");
h_reco_fakecorr->Draw("same");
leg3->Draw("same");
TH1D * hreco_fakecorr_rat = (TH1D*)h_reco_fakecorr->Clone("hreco_fakecorr_rat");
hreco_fakecorr_rat->Divide(h_reco_matched);

c3->cd(2);
c3->cd(2)->SetLogx();
drawClosure(l,hreco_fakecorr_rat,1);

c3->SaveAs("compare_reco_fake_corr.png");
c3->SaveAs("compare_reco_fake_corr.pdf");



//********************************************************************
//full correction
TFile * outfile = new TFile("pt_closure_PbPb.root","recreate");

TH1D * h_reco_fakecorr_effcorr = new TH1D("h_reco_fakecorr_effcorr",";p_{T};Arbitrary Units",ny,x);
nt_track->Draw("pt>>h_reco_fakecorr_effcorr","((1-fake)/eff)*weight*(trackselect && pt>0.5)"); 

h_reco_fakecorr_effcorr->SetMarkerColor(kRed);
h_reco_fakecorr_effcorr->SetLineColor(kRed);
TLegend *leg4 = new TLegend(0.15,0.1,0.8,0.4);
legFormat(leg4);
leg4->AddEntry(h_gen,"gen","p");
leg4->AddEntry(h_reco,"reco","p");
leg4->AddEntry(h_reco_fakecorr_effcorr,"reco corr","p");

TCanvas *c4 = new TCanvas("c4","",600,600);
makeMultiPanelCanvas(c4,1,2,0.0,0.0,0.15,0.15,0.02);

c4->cd(1);
c4->cd(1)->SetLogx();
c4->cd(1)->SetLogy();
h_gen->Draw();
h_gen->Write();
h_reco->Draw("same");
h_reco->Write();
TH1D * reco_over_gen = (TH1D*)h_reco->Clone("reco_over_gen");
reco_over_gen->Divide(h_gen);
reco_over_gen->Write();
h_reco_fakecorr_effcorr->Draw("same");
h_reco_fakecorr_effcorr->Write();
leg4->Draw("same");
TH1D * h_genreco_fullcorr=(TH1D*)h_reco_fakecorr_effcorr->Clone("h_genreco_fullcorr");
h_genreco_fullcorr->Divide(h_gen);
h_genreco_fullcorr->Write();

c4->cd(2);
c4->cd(2)->SetLogx();
drawClosure(l,h_genreco_fullcorr,2);

c4->SaveAs("compare_select_fullcorr.png");
c4->SaveAs("compare_select_fullcorr.pdf");
//outfile->Close();
}
