#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "TStyle.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include <sstream>
#include <vector>
const int nVAR = 6;

int StandardSatPlots(){
  gStyle->SetTitleW(0.95);
  gStyle->SetPadRightMargin(0.055);
  gStyle->SetPadLeftMargin(0.14);
  TCanvas *c = new TCanvas("c","c",0,0,1900,900);
  c->Divide(3,2);
  double size = 0.2;
  c->cd(1);
  //const char *file = "PEMetalonBlocked20degRet100degHeaterMirrorNofieldMirror.txt";
  const char *file = "Diode1stMeasOvernightPEM90degAn2Crossed072021.txt";//OvernightPEMetalonBlocked20degRet100degHeaterMirrorNofieldMirror.txt";
  //const char *file = "RTD105_4T_SilveredMirror_TempGauge05172021";
  //const char *file ="48hrs4.15T_PEM45degReflectedArmRet160deg051521.txt";
  //const char* file = "LargeTempChange1hrField4.15TPEMreflectedArmRet105A1_0degA2_88.9degHeaterOn051321.txt";
  //const char* file ="Overnight4.15TPEMreflectedArmRet80A1_0degA2_88.9degHeaterOnOff051321.txt";
  //const char* file = "Field4.15TPEMreflectedArmRet105A1_0degA2_88.9degHeaterOnOff051221.txt";//"OvernightRamp2to4TthenstableCameraPickoff050521.txt";
    //"OvernightRampto4.15TthenSteadyRet105_043021.txt";//"48hrscamerainfrontoflaserPEMon105deg042521.txt";//"Ramp5Tto1.076T04062021.txt";
  //const char* file = "RampDown4to0TPEMoutofbeamRet105deg042021.txt";
  //const char* file = "Overnight4T_Ret105TableThermOnWoodsupportBox042021.txt";
  //const char* file = "Overnight4T_Ret105TableThermOnPEMelectronicBox041921.txt";
  //const char* file = "Overnight4T_Ret105TableThermOnPEMelectronicBox041721.txt";
  //"15minRun3CoverMoves4.15TwithPEMoutofbeamRunningwithRet180degA284degTableThermOnPEMelectronicBox041921.txt";
  //"Weekend1.55TTempHallProbeLIA032921.txt";
  TGraph *gr = new TGraph((char*)file,"%lg, %*lg, %lg,");
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(size);
  gr->SetMarkerColor(kBlue);
  gr->Draw("ap");
  gPad->Update();
  gr->SetTitle("Measured 1F vs Time  (GF Fe foil 99.99%)");
  gr->GetXaxis()->SetTitle("Time (s)");
  gr->GetYaxis()->SetTitle("1F (V)");
  c->cd(2);
  TGraph *gr1 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %lg");
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(size);
  gr1->SetMarkerColor(kBlue);
  gr1->Draw("ap");
  gPad->Update();
  gr1->SetTitle("Measured Aux1 vs Time (GF Fe foil 99.99%)");
  gr1->GetXaxis()->SetTitle("Time (s)");
  gr1->GetYaxis()->SetTitle("Aux1 (DC Signal) (V)");
  c->cd(4);
  TGraph *gr2 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(size);
  gr2->SetMarkerColor(kBlue);
  gr2->Draw("ap");
  gPad->Update();
  gr2->SetTitle("Measured 2F vs Time (GF Fe foil 99.99%)");
  gr2->GetXaxis()->SetTitle("Time (s)");
  gr2->GetYaxis()->SetTitle("2F (V)");
  c->cd(5);
  TGraph *gr3 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerSize(size);
  gr3->SetMarkerColor(kBlue);
  gr3->Draw("ap");
  gPad->Update();
  gr3->SetTitle("Measured 1F/2F vs Time (GF Fe foil 99.99%)");
  gr3->GetXaxis()->SetTitle("Time (s)");
  gr3->GetYaxis()->SetTitle("1F/2F");
  c->cd(6);
  TGraph *gr5 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr5->SetMarkerStyle(8);
  gr5->SetMarkerSize(size);
  gr5->SetMarkerColor(kBlue);
  gr5->Draw("ap");
  gPad->Update();
  gr5->SetTitle("Measured B-field vs Time");
  gr5->GetXaxis()->SetTitle("Time (s)");
  gr5->GetYaxis()->SetTitle("B-field (T)");
  c->cd(3);
  TMultiGraph *mg = new TMultiGraph();
  TGraph *gr4 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr4->SetMarkerStyle(8);
  gr4->SetMarkerSize(size);
  gr4->SetMarkerColor(kBlue);
  gr4->Draw("ap");
  gPad->Update();
  
  mg->SetTitle("Measured Temperature (Red=Table, Blue=Hall) vs Time");
  mg->GetXaxis()->SetTitle("Time (s)");
  mg->GetYaxis()->SetTitle("Temperature (degC)");

  TGraph *gr6 = new TGraph((char*)file,"%lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr6->SetMarkerStyle(8);
  gr6->SetMarkerSize(size);
  gr6->SetMarkerColor(kRed);
  mg->Add(gr4);
  mg->Add(gr6);
  mg->Draw("ap");
  gPad->Update();
  
  TCanvas *c2 = new TCanvas("c2","c2",0,0,1200,900);
  c2->Divide(2,2);
  c2->cd(1);
  TGraph *gr7 = new TGraph((char*)file,"%*lg, %*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg,%*lg,%*lg, %lg");
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerSize(size);
  gr7->SetMarkerColor(kBlue);
  for(int i=0;i<gr7->GetN();++i){
    double x,y;
    gr7->GetPoint(i,x,y);
    gr7->SetPoint(i,y,x);
  }
  gr7->Draw("ap");
  gPad->Update();
  gr7->SetTitle("1F vs Table Temperature");
  gr7->GetXaxis()->SetTitle("Table Temp (degC)");
  gr7->GetYaxis()->SetTitle("1F (V)");
  c2->cd(2);
  TGraph *gr8 = new TGraph((char*)file,"%*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr8->SetMarkerStyle(8);
  gr8->SetMarkerSize(size);
  gr8->SetMarkerColor(kRed);
  for(int i=0;i<gr8->GetN();++i){
    double x,y;
    gr8->GetPoint(i,x,y);
    gr8->SetPoint(i,y,x);
  }
  gr8->Draw("ap");
  gPad->Update();
  gr8->SetTitle("1F/2F Ratio vs Table Temperature");
  gr8->GetXaxis()->SetTitle("Table Temp (degC)");
  gr8->GetYaxis()->SetTitle("1F/2F");
  c2->cd(3);
  TGraph *gr7h = new TGraph((char*)file,"%*lg, %*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg,%*lg,%lg, %*lg");
  gr7h->SetMarkerStyle(8);
  gr7h->SetMarkerSize(size);
  gr7h->SetMarkerColor(kBlue);
  for(int i=0;i<gr7h->GetN();++i){
    double x,y;
    gr7h->GetPoint(i,x,y);
    gr7h->SetPoint(i,y,x);
  }
  gr7h->Draw("ap");
  gPad->Update();
  gr7h->SetTitle("1F vs Hall Temperature");
  gr7h->GetXaxis()->SetTitle("Hall Temp (degC)");
  gr7h->GetYaxis()->SetTitle("1F (V)");
  c2->cd(4);
  TGraph *gr8h = new TGraph((char*)file,"%*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg, %*lg, %*lg, %*lg, %lg");
  gr8h->SetMarkerStyle(8);
  gr8h->SetMarkerSize(size);
  gr8h->SetMarkerColor(kRed);
  for(int i=0;i<gr8h->GetN();++i){
    double x,y;
    gr8h->GetPoint(i,x,y);
    gr8h->SetPoint(i,y,x);
  }
  gr8h->Draw("ap");
  gPad->Update();
  gr8h->SetTitle("1F/2F Ratio vs Hall Temperature");
  gr8h->GetXaxis()->SetTitle("Hall Temp (degC)");
  gr8h->GetYaxis()->SetTitle("1F/2F");
   
  TCanvas *c3 = new TCanvas("c3","c3",0,0,1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  TGraph *gr9 = new TGraph((char*)file,"%*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr9->SetMarkerStyle(8);
  gr9->SetMarkerSize(size);
  gr9->SetMarkerColor(kBlue);
  for(int i=0;i<gr9->GetN();++i){
    double x,y;
    gr9->GetPoint(i,x,y);
    gr9->SetPoint(i,y,x*1000);
  }
  gr9->Draw("ap");
  gPad->Update();
  gr9->SetTitle("2F vs Hall Temperature");
  gr9->GetXaxis()->SetTitle("Hall Temp (degC)");
  gr9->GetYaxis()->SetTitle("2F (mV)");
  gr9->Fit("pol1");
  c3->cd(2);
  TGraph *gr10 = new TGraph((char*)file,"%*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr10->SetMarkerStyle(8);
  gr10->SetMarkerSize(size);
  gr10->SetMarkerColor(kRed);
  for(int i=0;i<gr10->GetN();++i){
    double x,y;
    gr10->GetPoint(i,x,y);
    gr10->SetPoint(i,y,x*1000);
  }
  gr10->Draw("ap");
  gPad->Update();
  gr10->SetTitle("2F vs Table Temperature");
  gr10->GetXaxis()->SetTitle("Table Temp (degC)");
  gr10->GetYaxis()->SetTitle("2F (mV)");
  cout<<gr9->GetCorrelationFactor()<<" "<<gr10->GetCorrelationFactor()<<endl;
  gStyle->SetOptFit(1111);
  gr10->Fit("pol1");
  return 0;
}
