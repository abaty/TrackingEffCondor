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
#include "../condor_trk_corr/ntupler/trackTree.C"

void reweight(){
 TH1D::SetDefaultSumw2();

 //converted to nb
 // full stats 

  float crossSection[6] = {1.075E-02,1.025E-03,9.865E-05,3.069E-05,1.129E-05,0}; 

  const int nbins =   5;
  float bins[nbins+1]= {30,50,80,100,120,1000};
  float weight[nbins]={};

  TCanvas * c2 = new TCanvas("c2","c2",700,500);
  c2->SetLogy();

  TH1F * pthatDist = new TH1F("pthatDist","pthatDist",nbins,bins);
  TH1F * pthatDistFine = new TH1F("pthatDistFine","pthatDistFine",100,0,500);
  TH1F * vzDistMC = new TH1F("vzDistMC","vzDistMC",200,-50,50);

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
   //ftrk[ifile] = new trackTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile]));
   fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile]));
   evtSelFile[ifile] = new TFile(Form("%s/%s.root",directory.Data(),infname[ifile]),"read");
   evtSel[ifile] = (TTree*) evtSelFile[ifile]->Get("skimanalysis/HltTree");
   evtSel[ifile]->SetBranchAddress("pcollisionEventSelection", &pcoll[ifile]); 
  }

  int nevents = 150000;
//file loop here
  for(int ifile=2; ifile<5; ifile++){
  std::cout<<ifile<<std::endl;
  int nentries = fhi[ifile]->GetEntriesFast();

//event loop 
  if(nevents<nentries) nentries = nevents;
  for(int jentry=0;jentry<nentries;jentry++){
  if((jentry%10000)==0) std::cout<<jentry<<"/"<<nentries<< "   File:" << ifile <<std::endl;

  //ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

//vertexShift can be used to move the center of the vertex distribution if needed
  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;
 
  pthatDist->Fill(fjet[ifile]->pthat);
  pthatDistFine->Fill(fjet[ifile]->pthat);
  }
 }

pthatDist->Draw("h");
c2->SaveAs("pthatDist.png");
pthatDistFine->Draw("h");
c2->SaveAs("pthatDistFine.png");

TH1F * weights = new TH1F("weights","weights",nbins,bins);
  using namespace std;
  for(int i=2; i<nbins; i++){
    weight[i]=(crossSection[i]-crossSection[i+1])*1000000/((float)pthatDist->GetBinContent(i+1));
    weights->SetBinContent(i+1,weight[i]);

    cout << bins[i] << "-" << bins[i+1] << " GeV/c"<< endl;
    cout << "difference in cross section: " << crossSection[i]-crossSection[i+1] <<endl;
    cout << "total number of events:      " << pthatDist->GetBinContent(i+1) << endl;
    cout << "weight:                      " << weight[i] << "\n"<< endl;
  }

  cout << "{" << weight[0] <<","<<weight[1]<<","<<weight[2]<<","<<weight[3]<<","<<weight[4]<< "}" << endl;

  pthatDist->Multiply(weights);
  pthatDist->Draw("h");
  c2->SaveAs("pthatDistWeighted.png");

  TH1F * weightsFine = new TH1F("weightsFine","weightsFine",100,0,500);
  for(int i=0;i<100;i++){
    //if(5*(i) <50)        weightsFine->SetBinContent(i+1,weight[0]);
    //else if(5*(i)<80)  weightsFine->SetBinContent(i+1,weight[1]);
    //else if(5*(i)<100) weightsFine->SetBinContent(i+1,weight[2]);
    if(5*(i)<100) weightsFine->SetBinContent(i+1,weight[2]);
    else if(5*(i)<120) weightsFine->SetBinContent(i+1,weight[3]);
    else                 weightsFine->SetBinContent(i+1,weight[4]);
  }

  pthatDistFine->Multiply(weightsFine);
  pthatDistFine->Draw("h");
  c2->SaveAs("pthatDistWeightedFine.png");
//





///***********************************************************cent reweighting part******************************************
TCanvas * c3 = new TCanvas("c3","c3",500,700);

for(int ifile=2; ifile<5; ifile++){
  std::cout<<ifile<<std::endl;
  int nentries = fhi[ifile]->GetEntriesFast();

//event loop
  if(nevents<nentries) nentries = nevents;
  for(int jentry=0;jentry<nentries;jentry++){
  if((jentry%10000)==0) std::cout<<jentry<<"/"<<nentries<< "   File:" << ifile <<std::endl;

  //ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

//vertexShift can be used to move the center of the vertex distribution if needed
  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;

  //if(fjet[ifile]->pthat<50) vzDistMC->Fill(vz,weight[0]);
  //else if(fjet[ifile]->pthat<80) vzDistMC->Fill(vz,weight[1]);
  //else if(fjet[ifile]->pthat<100) vzDistMC->Fill(vz,weight[2]);
  if(fjet[ifile]->pthat<100) vzDistMC->Fill(vz,weight[2]);
  else if(fjet[ifile]->pthat<120) vzDistMC->Fill(vz,weight[3]);
  else  vzDistMC->Fill(vz,weight[4]);
  }
 }

  TFile * dataFile = new TFile("/mnt/hadoop/cms/store/user/velicanu/HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged/0.root","read");
  TTree * data = (TTree*) dataFile->Get("hiEvtAnalyzer/HiTree");
  TTree * dataFriend =(TTree*) dataFile->Get("skimanalysis/HltTree");
  data->AddFriend(dataFriend); 
 
  TH1F * dataHist = new TH1F("dataHist","dataHist",100,0,200);
  data->Draw("hiBin>>dataHist","abs(vz)<15 && pcollisionEventSelection");
  dataHist->Scale(1/(dataHist->Integral(1,100,"width")));
  dataHist->Draw("h");
  c3->SaveAs("data_cent_dist.png");

  TH1F * vzHistdata = new TH1F("vzHistdata","vzHistdata",200,-50,50);
  data->Draw("vz>>vzHistdata","");

  float vertexShift = vzDistMC->GetMean()-vzHistdata->GetMean();

  std::cout << "<vz>MC: " << vzDistMC->GetMean() << std::endl;
  std::cout << "<vz>data: " << vzHistdata->GetMean() << std::endl;
  std::cout << "vertexShift: " << vertexShift << std::endl;

TH1F * MC = new TH1F("MC","MC",100,0,200);


 for(int ifile=2; ifile<5; ifile++){
  std::cout<<ifile<<std::endl;
  int nentries = fhi[ifile]->GetEntriesFast();

//event loop
  if(nevents<nentries) nentries = nevents;
  for(int jentry=0;jentry<nentries;jentry++){
  if((jentry%10000)==0) std::cout<<jentry<<"/"<<nentries<< "   File:" << ifile <<std::endl;

  //ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  evtSel[ifile]->GetEntry(jentry);

//vertexShift can be used to move the center of the vertex distribution if needed
  float cent=fhi[ifile]->hiBin;
  float vz = fhi[ifile]->vz;

  if(abs(vz-vertexShift)>15 || !(pcoll[ifile])) continue;
  if(fjet[ifile]->pthat<50) MC->Fill(cent,weight[0]);
  else if(fjet[ifile]->pthat<80) MC->Fill(cent,weight[1]);
  else if(fjet[ifile]->pthat<100) MC->Fill(cent,weight[2]);
  else if(fjet[ifile]->pthat<120) MC->Fill(cent,weight[3]);
  else  MC->Fill(cent,weight[4]);
  }
 }
  
  double integral = MC->Integral(1,100,"width");
  MC->Scale(1/(integral));
  MC->Draw("h");
  c3->SaveAs("MC_cent_dist.png");

  dataHist->Divide(MC);
  dataHist->Draw("h");
  c3->SaveAs("centrality_weight.png");

  TFile * outf = new TFile("centrality_weights.root","recreate");
  dataHist->SetName("centrality_weight");
  dataHist->SetTitle("centrality_weight");
  dataHist->Write();
  outf->Close();

  MC->Scale(integral);
  MC->Multiply(dataHist);
  MC->Draw("h");
  c3->SaveAs("MC_cent_dist_weighted.png");


}
