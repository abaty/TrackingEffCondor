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
  histo->GetYaxis()->SetLabelSize(13);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleSize(20);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetTitleOffset(1.5);
  histo->GetYaxis()->CenterTitle();
  histo->GetXaxis()->SetLabelSize(13);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetXaxis()->SetTitleSize(20);
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleOffset(1.5);
}

void legFormat(TLegend* leg, float low, float high){
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0,"CMS Preliminary","");
  leg->AddEntry((TObject*)0,Form("PYTHIA+HYDJET, %f<p_{T}<%f",low,high),"");
  leg->AddEntry((TObject*)0,"VsCalo jets, HI tracking","");
}

void drawClosure(TLine * l, TH1D * histo){
histo->SetMaximum(1.5);
histo->SetMinimum(0.5);
histo->GetYaxis()->SetTitle("Closure");
histo->Draw();
fixAxes(histo);
l->Draw("same");
}

void pt_closure_test(int m = 0){

if(m<0 || m>3){
  m=0;
  std::cout << "invalid mode! Assuming cent mode." << std::endl;
}

const char* var[4] = {"cent","rmin","eta","phi"};
const char* label[4] = {"cent(%)","r_{min}","#eta","#phi"};
const int bins[4] = {50,50,50,50};
const float lowerBin[4] = {0,0,-2.4,-TMath::Pi()};
const float higherBin[4] = {100,5,2.4,TMath::Pi()};
const float x1[4]= {0.55  ,0.3, 0.5 ,0.55};
const float y1[4]= {0.65 ,0.1, 0.65 ,0.65};
const float x2[4]= {0.95 ,0.8, 0.9 ,0.95};
const float y2[4]= {0.95 ,0.4, 0.95 ,0.95};
float pt_low = 0.85;
float pt_high= 1.0; 


TH1D::SetDefaultSumw2();
TFile * f= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_pthatCombo_100ksmall.root");
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
double delta=(log(bin_pt_max)-log(bin_pt_min))/(50*log(10));
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix); 
 inix+=delta;
}

//eff correction
TH1D * h_gen = new  TH1D("h_gen",Form(";%s;N",label[m]),bins[m],lowerBin[m],higherBin[m]);
TH1D * h_gen_select = new  TH1D("h_gen_select",Form(";%s;N_{evt}",label[m]),bins[m],lowerBin[m],higherBin[m]);
TH1D * h_gen_matched_select_corr = new TH1D("h_gen_matched_select_corr",Form(";%s;N_{evt}",label[m]),bins[m],lowerBin[m],higherBin[m]);

nt_particle->Draw(Form("%s>>h_gen",var[m]),Form("weight*(pt>%f && pt<%f)",pt_low,pt_high));
nt_particle->Draw(Form("%s>>h_gen_select",var[m]),Form("weight*(trackselect && pt>%f && pt<%f)",pt_low,pt_high));
nt_particle->Draw(Form("%s>>h_gen_matched_select_corr",var[m]),Form("(1/eff)*weight*(trackselect && pt>%f && pt<%f)",pt_low,pt_high));

h_gen->SetMarkerColor(1);
h_gen->SetMarkerStyle(25);
h_gen->SetLineWidth(1);
h_gen_select->SetMarkerColor(1);
h_gen_matched_select_corr->SetMarkerColor(kRed);
h_gen_matched_select_corr->SetLineColor(kRed);

TLegend *leg2;
if(m==0) leg2 = new TLegend(x1[m],y1[m],x2[m],y2[m]); 
if(m==1) leg2 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==2) leg2 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==3) leg2 = new TLegend(x1[m],y1[m],x2[m],y2[m]);

legFormat(leg2,pt_low,pt_high);
leg2->AddEntry(h_gen,"gen","p");
leg2->AddEntry(h_gen_select,"matched gen","p");
leg2->AddEntry(h_gen_matched_select_corr,"corrected matched gen","p");

TCanvas * c2 = new TCanvas("c2","",600,600);
makeMultiPanelCanvas(c2,1,2,0.0,0.0,0.15,0.15,0.02);
c2->cd(1);
if(m<2) c2->cd(1)->SetLogy();
else{
  h_gen->SetMaximum(1.6*h_gen->GetBinContent(h_gen->GetMaximumBin()));
  h_gen->SetMinimum(0.7*h_gen_select->GetBinContent(h_gen_select->GetMinimumBin()));
}

h_gen->Draw();
fixAxes(h_gen);
h_gen_select->Draw("same");
h_gen_matched_select_corr->Draw("same");
leg2->Draw("same");

TH1D * hgen_corr_rat = (TH1D*)h_gen_matched_select_corr->Clone("hgen_corr_rat");
hgen_corr_rat->Divide(h_gen);
c2->cd(2);

TLine * l; 
if(m==0) l = new TLine(lowerBin[m],1,higherBin[m],1);
if(m==1) l = new TLine(lowerBin[m],1,higherBin[m],1);
if(m==2) l = new TLine(lowerBin[m],1,higherBin[m],1);
if(m==3) l = new TLine(lowerBin[m],1,higherBin[m],1);

drawClosure(l,hgen_corr_rat);

c2->SaveAs(Form("compare_gen_select_corr_%spt%d%d.png",var[m],(int)(10*pt_low),(int)(10*pt_high)));
c2->SaveAs(Form("compare_gen_select_corr_%spt%d%d.pdf",var[m],(int)(10*pt_low),(int)(10*pt_high)));

//****************************************************************************
//fake correction
TH1D * h_reco = new  TH1D("h_reco",Form(";%s;N",label[m]),bins[m],lowerBin[m],higherBin[m]);
TH1D * h_reco_matched = new  TH1D("h_reco_matched",Form(";%s;N_{evt}",label[m]),bins[m],lowerBin[m],higherBin[m]);
TH1D * h_reco_fakecorr = new TH1D("h_reco_fakecorr",Form(";%s;N_{evt}",label[m]),bins[m],lowerBin[m],higherBin[m]);
nt_track->Draw(Form("%s>>h_reco",var[m]),Form("weight*(trackselect && pt>%f && pt<%f)",pt_low,pt_high)); 
nt_track->Draw(Form("%s>>h_reco_matched",var[m]),Form("weight*(trackselect && !trkfake && pt>%f && pt<%f)",pt_low,pt_high)); 
nt_track->Draw(Form("%s>>h_reco_fakecorr",var[m]),Form("(1-fake)*weight*(trackselect && pt>%f && pt<%f)",pt_low,pt_high)); 

h_reco_matched->SetMarkerColor(1);
h_reco_matched->SetMarkerStyle(25);
h_reco_fakecorr->SetMarkerColor(kRed);
h_reco_fakecorr->SetLineColor(kRed);

TLegend *leg3;
if(m==0) leg3 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==1) leg3 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==2) leg3 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==3) leg3 = new TLegend(x1[m],y1[m],x2[m],y2[m]);

legFormat(leg3,pt_low,pt_high);
leg3->AddEntry(h_reco,"reco","p");
leg3->AddEntry(h_reco_matched,"matched reco","p");
leg3->AddEntry(h_reco_fakecorr,"fake corrected reco","p");

TCanvas * c3 = new TCanvas("c3","",600,600);
makeMultiPanelCanvas(c3,1,2,0.0,0.0,0.15,0.15,0.02);
c3->cd(1);
if(m<2) c3->cd(1)->SetLogy();
else{
  h_reco->SetMaximum(1.6*h_reco->GetBinContent(h_reco->GetMaximumBin()));
  h_reco->SetMinimum(0.7*h_reco_matched->GetBinContent(h_reco_matched->GetMinimumBin()));
}

h_reco->Draw();
fixAxes(h_reco);

h_reco_matched->Draw("same");
h_reco_fakecorr->Draw("same");
leg3->Draw("same");
TH1D * hreco_fakecorr_rat = (TH1D*)h_reco_fakecorr->Clone("hreco_fakecorr_rat");
hreco_fakecorr_rat->Divide(h_reco_matched);

c3->cd(2);
drawClosure(l,hreco_fakecorr_rat);

c3->SaveAs(Form("compare_reco_fake_corr_%spt%d%d.png",var[m],(int)(10*pt_low),(int)(10*pt_high)));
c3->SaveAs(Form("compare_reco_fake_corr_%spt%d%d.pdf",var[m],(int)(10*pt_low),(int)(10*pt_high)));



//********************************************************************
//full correction
TH1D * h_reco_fakecorr_effcorr = new TH1D("h_reco_fakecorr_effcorr",Form(";%s;N",label[m]),bins[m],lowerBin[m],higherBin[m]);
nt_track->Draw(Form("%s>>h_reco_fakecorr_effcorr",var[m]),Form("((1-fake)/eff)*weight*(trackselect && pt>%f && pt<%f)",pt_low,pt_high)); 

h_reco_fakecorr_effcorr->SetMarkerColor(kRed);
h_reco_fakecorr_effcorr->SetLineColor(kRed);
TLegend *leg4;
if(m==0) leg4 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==1) leg4 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==2) leg4 = new TLegend(x1[m],y1[m],x2[m],y2[m]);
if(m==3) leg4 = new TLegend(x1[m],y1[m],x2[m],y2[m]);

legFormat(leg4,pt_low,pt_high);
leg4->AddEntry(h_gen,"gen","p");
leg4->AddEntry(h_reco,"reco","p");
leg4->AddEntry(h_reco_fakecorr_effcorr,"reco corr","p");

TCanvas *c4 = new TCanvas("c4","",600,600);
makeMultiPanelCanvas(c4,1,2,0.0,0.0,0.15,0.15,0.02);

c4->cd(1);
if(m<2) c4->cd(1)->SetLogy();
h_gen->Draw();
h_reco->Draw("same");
h_reco_fakecorr_effcorr->Draw("same");
leg4->Draw("same");
TH1D * h_genreco_fullcorr=(TH1D*)h_reco_fakecorr_effcorr->Clone("h_genreco_fullcorr");
h_genreco_fullcorr->Divide(h_gen);

c4->cd(2);
drawClosure(l,h_genreco_fullcorr);

c4->SaveAs(Form("compare_select_fullcorr_%spt%d%d.png",var[m],(int)(10*pt_low),(int)(10*pt_high)));
c4->SaveAs(Form("compare_select_fullcorr_%spt%d%d.pdf",var[m],(int)(10*pt_low),(int)(10*pt_high)));
}
