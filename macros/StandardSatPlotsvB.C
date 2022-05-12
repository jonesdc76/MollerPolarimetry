#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include <sstream>
#include <vector>
const int nVAR = 6;

int StandardSatPlotsvB(){
  gStyle->SetTitleW(0.95);
  gStyle->SetPadRightMargin(0.055);
  gStyle->SetPadLeftMargin(0.14);
  TCanvas *c = new TCanvas("c","c",0,0,1900,900);
  c->Divide(3,2);
  const char *filename = "Rampto4TPEMreflectedArmRet105A10degA288.9deg051221.txt";//"Ramp5to2TPickoffToCamera050421.txt";
  //const char *filename = "RampDown2ndGF99.99033021.txt";
  //const char *filename = "RampUpto4.8T2nfGF99.99An2_134.5deg040521.txt";
  //const char *filename = "RampDown4.8to2.1TA2_134.5deg040521.txt";
  //const char *filename = "RampUp2.8to5TGF3_An2_105deg_040621.txt";
  //const char *filename = "GF1Repol_Ramp1to5TAn2_162deg_Ret110040821.txt";
  //const char *filename = "RampDown4to0TPEMoutofbeamRet105deg042021.txt";//"Ramp1.038to4.4115FeFoilInSlot3_03302021.txt";
  //const char *filename = "RampUp2.8to5TGF3_An2_105deg_040621.txt";
  //const char *filename = "Ramp5Tto1.076T04062021.txt";
  c->cd(1);
  TGraph *gr = new TGraph(filename,"%*lg, %lg, %lg,");
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.1);
  gr->SetMarkerColor(kBlue);
  gr->Draw("ap");
  gPad->Update();
  gr->SetTitle("Measured 1F vs B-field (GF Fe foil 99.99%)");
  gr->GetXaxis()->SetTitle("B-field (T)");
  gr->GetYaxis()->SetTitle("1F (V)");
  c->cd(2);
  TGraph *gr1 = new TGraph(filename,"%*lg, %lg, %*lg, %*lg, %lg");
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(0.1);
  gr1->SetMarkerColor(kBlue);
  gr1->Draw("ap");
  gPad->Update();
  gr1->SetTitle("Measured Aux1 vs B-field (GF Fe foil 99.99%)");
  gr1->GetXaxis()->SetTitle("B-field (T)");
  gr1->GetYaxis()->SetTitle("Aux1 (DC Signal) (V)");
  c->cd(4);
  TGraph *gr2 = new TGraph((char*)filename,"%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(0.1);
  gr2->SetMarkerColor(kBlue);
  gr2->Draw("ap");
  gPad->Update();
  gr2->SetTitle("Measured 2F vs B-field (GF Fe foil 99.99%)");
  gr2->GetXaxis()->SetTitle("B-field (T)");
  gr2->GetYaxis()->SetTitle("2F (V)");
  c->cd(5);
  TGraph *gr3 = new TGraph((char*)filename,"%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerSize(0.1);
  gr3->SetMarkerColor(kBlue);
  gr3->Draw("ap");
  gPad->Update();
  gr3->SetTitle("Measured 1F/2F vs B-field (GF Fe foil 99.99%)");
  gr3->GetXaxis()->SetTitle("B-field (T)");
  gr3->GetYaxis()->SetTitle("1F/2F");
  c->cd(6);
  TGraph *gr5 = (TGraph*)gr3->Clone("gr5");
  gr5->Draw("ap");
  gr5->GetXaxis()->SetRangeUser(2.7,6);
  //gr5->GetYaxis()->SetRangeUser(0.02,0.0216);
  // TGraph *gr5 = new TGraph((char*)filename,"%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  // gr5->SetMarkerStyle(8);
  // gr5->SetMarkerSize(0.1);
  // gr5->SetMarkerColor(kBlue);
  // gr5->Draw("ap");
  // gPad->Update();
  // gr5->SetTitle("Measured Fringe B-field vs Fringe B-field");
  // gr5->GetXaxis()->SetTitle("B-field (kG)");
  // gr5->GetYaxis()->SetTitle("B-field (kG)");
  c->cd(3);
  TMultiGraph *grm = new TMultiGraph();
  grm->SetTitle("Measured Temperature (Hall=Blue, Optics table = Red)  vs B-field");
  TGraph *gr4 = new TGraph((char*)filename,"%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr4->SetMarkerStyle(8);
  gr4->SetMarkerSize(0.1);
  gr4->SetMarkerColor(kBlue);
  grm->Add(gr4);

  TGraph *gr6 = new TGraph((char*)filename,"%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr6->SetMarkerStyle(8);
  gr6->SetMarkerSize(0.1);
  gr6->SetMarkerColor(kRed);
  grm->Add(gr6);
  grm->Draw("ap");
  gPad->Update();
  grm->GetXaxis()->SetTitle("B-field (T)");
  grm->GetYaxis()->SetTitle("Temperature (degC)");
  grm->Draw("ap");
  gPad->Update();
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,450);
  gr3->Draw("ap");

  return 0;
}
