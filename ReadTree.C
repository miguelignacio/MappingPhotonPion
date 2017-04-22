#include "TFile.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TRandom.h"
using namespace std;

double Resolution(double E){
	double sigma = E*TMath::Sqrt( pow(0.081/TMath::Sqrt(E),2) + pow(0.021,2));
	return sigma; 
}

void ReadTree(){

  TRandom *RANDOM = new TRandom();
  std::cout<<"Starting Analysis" << std::endl;
  auto myFile =  TFile::Open("Pi0Events_v2.root","READ");

  if (!myFile || myFile->IsZombie()) {
    return;
  }
  myFile->Print();

  auto tout =(TTree*)myFile->Get("T");
  if(!tout) {
    std::cout <<"Problem getting TTree" << std::endl;
    return;
  }
  
  TTreeReader myReader(tout);
  TTreeReaderValue<std::vector<int>> ID(myReader, "ID");
  TTreeReaderValue<std::vector<float>> Pt(myReader, "Pt");
  TTreeReaderValue<std::vector<float>> Eta(myReader, "Eta");
  TTreeReaderValue<std::vector<float>> Phi(myReader, "Phi");
  TTreeReaderValue<std::vector<float>> Theta(myReader, "Theta");
  TTreeReaderValue<std::vector<float>> E(myReader, "E");
  
  int nevent = 0;
  
  int nbins_pTPion = 200;
  int nbins_pTPhoton = 40;
  int nbins_MassPion = 100;
  
  double min_pTPion = 0.0;
  double max_pTPion= 20.0;
  double min_pTPhoton = 0.0;
  double max_pTPhoton = 20.0;
  
  
  
  
  double binwidth = max_pTPhoton/nbins_pTPhoton;
  
  Int_t bins[4]    = {nbins_pTPion, nbins_pTPhoton, nbins_pTPhoton, nbins_MassPion};
  Double_t xmin[4] = {min_pTPion,    min_pTPhoton,  min_pTPhoton, 0.0};
  Double_t xmax[4] = {max_pTPion,    max_pTPhoton,  max_pTPhoton, 0.200};

  THnSparseD* h = new THnSparseD("h", "#pT pion; # pT Photon; pT Photon; Mass ",4,bins,xmin,xmax);
  h->Sumw2();
  while (myReader.Next()) {
    nevent +=1;
    //if(nevent%1000==0) std::cout << " Event # " << nevent << std::endl;
 
   TLorentzVector photon_1, photon_2, photon_smeared_1, photon_smeared_2, vpion, vpion_smeared;
   photon_1.SetPtEtaPhiE(E->at(2)*TMath::Sin(Theta->at(2)), Eta->at(2), Phi->at(2), E->at(2));
   photon_2.SetPtEtaPhiE(E->at(3)*TMath::Sin(Theta->at(3)), Eta->at(3), Phi->at(3), E->at(3));
   
   //std::cout << E->at(2)*TMath::Cos(Theta->at(2)) << " " <<  E->at(2)*TMath::Sin(Theta->at(2)) << " " << Pt->at(2) << std::endl;
   //photon_1.SetPtEtaPhiE(Pt->at(2), Eta->at(2), Phi->at(2), E->at(2));
   //photon_2.SetPtEtaPhiE(Pt->at(3), Eta->at(3), Phi->at(3), E->at(3));
   vpion = photon_1 + photon_2;
   
   double dE_1 = RANDOM->Gaus(0, Resolution(photon_1.E()));
   double dE_2 = RANDOM->Gaus(0, Resolution(photon_2.E()));
   photon_smeared_1.SetPtEtaPhiE((E->at(2)+dE_1)*TMath::Sin(Theta->at(2)), Eta->at(2), Phi->at(2), E->at(2) + dE_1);
   photon_smeared_2.SetPtEtaPhiE((E->at(3)+dE_2)*TMath::Sin(Theta->at(3)), Eta->at(3), Phi->at(3), E->at(3) + dE_2);
   
   vpion_smeared = photon_smeared_1 + photon_smeared_2;
     
 
   double entries[4] = {Pt->at(1), photon_1.Pt(), photon_smeared_1.Pt(), vpion_smeared.M()};
   h->Fill(entries);
   double entries2[4] = {Pt->at(1), photon_2.Pt(), photon_smeared_2.Pt(), vpion_smeared.M()};
   h->Fill(entries2);
   } // end looping events
   
   
  TFile* fout = new TFile("fout.root","RECREATE");
  h->Write("Hsparse");
   
  
  }// end function
   
   