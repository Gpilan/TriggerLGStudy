#include "CBDsimRootInterface.h"
#include "CBDsimInterface.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TH3.h"
#include "TMath.h"
#include "TFile.h"


#include <iostream>
#include <string>
#include <cmath>


int main(int argc, char* argv[]){
  
  TString filename =argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  TH1F* tEdep = new TH1F("totEdep","Total Energy deposit;MeV;Evt",100,low,high);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);



 
  TH1F* tT = new TH1F("time","time;ns;p.e.",24000,0,240);
  tT->Sumw2(); tT->SetLineColor(kBlue); tT->SetLineWidth(2);

  TH1F* tWav = new TH1F("wavlen","wavelength;nm;p.e.",60,300.,900.);
  tWav->Sumw2(); tWav->SetLineColor(kRed); tWav->SetLineWidth(2);


  TH1F* tNhit = new TH1F("nHits","Number of p.e./SiPM;p.e.;n",200,0.,2000.);
  tNhit->Sumw2(); tNhit->SetLineColor(kRed); tNhit->SetLineWidth(2);

  TH1I* thits = new TH1I("p.e.","p.e./evt;p.e.;evt",100,0,1500);thits->Sumw2(); thits->SetLineColor(kBlack); thits->SetLineWidth(2);
  //std::cout<<"11"<<std::endl;

  TH1F* tTurnOnTime = new TH1F("risetime","risetime;t(ns);evt",24000,0,240);

  TH2I *tsipmXY = new TH2I("sipmXY","sipmXY;X;Y",30,0,30,30,0,30);

  TH2F *tXYriseTime = new TH2F("XYturnontime","turn on time by sipm;X;Y",30,0,30,30,0,30);

  TH1F* tenergy = new TH1F("energy detected", "energy detected;MeV;evt",60,low,high);

  TF1 *genergy = new TF1("genergy","gaus",10.3,10.7);genergy->SetLineColor(2);genergy->SetLineWidth(2);

  TH1I* thitstotal = new TH1I("p.e._total","p.e./evt(total);p.e.;evt",100,1000,3000);thitstotal->Sumw2(); thitstotal->SetLineColor(kBlack); thitstotal->SetLineWidth(2);

  TH2I* genparticle2D = new TH2I("2Dgenparticle","2D generated particle;x;y",30,0,30,30,0,30);
  TH2I* gengamma2D = new TH2I("2Dgengamma","generated gamma;x;y",30,0,30,30,0,30);
  TH2I* genelectron2D = new TH2I("2Dgenelectron","generated electron;x;y",30,0,30,30,0,30);
  TH2I* genpositron2D = new TH2I("2Dgenppositron","generated positron;x;y",30,0,30,30,0,30);
  TH1F* risetimesipm = new TH1F("","1sipm;t(ns);p.e.",24000,0,240);

  CBDsimRootInterface* cbdInterface = new CBDsimRootInterface(std::string(filename)+".root");
  cbdInterface->set();
  
  
  CBDsimInterface::CBDsimTimeStruct tmp;
  unsigned int entries =cbdInterface->entries();
  float totenergy=0.;
  int tothit=0;
  TCanvas* c1 = new TCanvas("c1","");
  int zedepnum=0;
  while(cbdInterface->numEvt()<entries){
    if (cbdInterface->numEvt()%50 ==0)printf("Analyzing %dth event ...\n", cbdInterface->numEvt());
    CBDsimInterface::CBDsimEventData cbdEvt;
    cbdInterface->read(cbdEvt);
    int nhit=0;
    float energy=0;

    float Edep =0;

    TH1F* turnontime= new TH1F("","1MeV 1evt;t(ns);p.e.",24000,0,240);

    for (auto edep : cbdEvt.Edeps){
      Edep += edep.Edep;
    }
    tEdep->Fill(Edep);
    {
      const CBDsimInterface::CBDsimTowerData* towers[] = {&cbdEvt.towerT1, &cbdEvt.towerT2};
      for (const auto* ptw : towers) {
        const auto& tower = *ptw;
      for (auto sipm:tower.SiPMs){
        tsipmXY->Fill(sipm.x,sipm.y,sipm.count);
        tNhit->Fill(sipm.count);
        int hit = sipm.count;
        nhit+=hit;
        for (const auto timepair:sipm.timeStruct){
          tT->Fill(timepair.first.first+0.005,timepair.second);
          turnontime->Fill(timepair.first.first+0.005,timepair.second);
          if(cbdInterface->numEvt()==100)risetimesipm->Fill(timepair.first.first+0.005,timepair.second);
        }
        for (const auto wavpair: sipm.wavlenSpectrum){
          tWav->Fill(wavpair.first.first,wavpair.second);
        }
        if(cbdInterface->numEvt()==100){
          int risetimesipmbin=risetimesipm->FindFirstBinAbove(0.5);
          if(Edep!=0) tXYriseTime->Fill(sipm.x,sipm.y,risetimesipm->GetBinCenter(risetimesipmbin));
          risetimesipm->Reset();
        }
      }
      }
    }
    for (auto physical:cbdEvt.totPhysicals){

     
      if (cbdInterface->numEvt()==100){
        double xphy =30-(physical.x-424.22)/5.052;
        double yphy = (physical.y+75.78)/5.052;
        double zphy = physical.z/10;
        genparticle2D->Fill(yphy,xphy);
        if(physical.particleID==0)gengamma2D->Fill(yphy,xphy,10);
        if(physical.particleID==1)genelectron2D->Fill(yphy,xphy,20);
        if(physical.particleID==2)genpositron2D->Fill(yphy,xphy,30);
        std::cout<<"particle ID is "<<physical.particleID<<"| particle position is "<<yphy<<","<<xphy<<"| particle energy is "<<physical.energy<<std::endl;

      }

      
    }
    tothit+=nhit;
    totenergy+=energy;
    thits->Fill(nhit);
    if(Edep!=0)thitstotal->Fill(nhit);
    if(Edep!=0) tenergy->Fill(nhit*0.000333838);//21.07.20 calib const 0.00033838
    //onecell 0.000333501 total 0.000322373
   
    int risetimebin=turnontime->FindFirstBinAbove(0.5);
    if(Edep!=0) tTurnOnTime->Fill(turnontime->GetBinCenter(risetimebin));

  }
  
  cbdInterface->close();

  TCanvas* c = new TCanvas("c","");
  tenergy->Fit(genergy,"R+");

  tXYriseTime->SetStats(kFALSE);
 
  thitstotal->Draw("Hist");c->SaveAs(filename+"_totalp.e._total.png");
  tTurnOnTime->Draw("Hist");c->SaveAs(filename+"_risetime.png");

  tenergy->Draw("Hist"); genergy->Draw("same"); c->SaveAs(filename+"_energy.png");
  tEdep->Draw("Hist"); c->SaveAs(filename+"_Edep.png");
  tT->Draw("Hist"); c->SaveAs(filename+"_t.png");
  tWav->Draw("Hist"); c->SaveAs(filename+"_wav.png");
  tNhit->Draw("Hist"); c->SaveAs(filename+"_hit.png");
  thits->Draw("Hist");  c->SaveAs(filename+"_totalp.e._back.png");
  genparticle2D->Draw("COL");c->SaveAs(filename+"_XY_gen.png");
  gengamma2D->Draw("COL");c->SaveAs(filename+"_XY_gamma.png");
  genelectron2D->Draw("COL");c->SaveAs(filename+"_XY_electron.png");
  genpositron2D->Draw("COL");c->SaveAs(filename+"_XY_positron.png");
  tXYriseTime->Draw("COLZ");c->SaveAs(filename+"_TOtimesipm.png");


  TFile* riseFile = new TFile(filename+"_risetime0.5.root","RECREATE");
  riseFile->WriteTObject(tTurnOnTime);
  riseFile->Close();
  
}
