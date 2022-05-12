#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include <iostream>
#include "TLegend.h"
#include "TPad.h"
#include "TAxis.h"

void norm(TGraph* gr, bool reverse = 0, double tweak = 1){
  double scale = 0;
  double n=5;
  int ni = 0;
  for(int i=reverse? 0:gr->GetN()-1;;--i){
    double x,y;
    gr->GetPoint(i,x,y);
    scale += y/n;
    ++ni;
    cout<<i<<" "<<ni<<" "<<y<<endl;
    if(ni>n-0.00001)break;
    if(reverse)i+=2;
  }
  for(int i=0;i<gr->GetN();++i){
    double x,y;
    gr->GetPoint(i,x,y);
    cout<<i<<" "<<x<<" "<<y<<endl;
    gr->SetPoint(i,x,y/scale*tweak);
  }
  cout<<scale<<endl;
}
void scale(TGraph* gr){
  double scale = 1000;
  double n=5;
  int ni = 0;
  for(int i=0;i<gr->GetN();++i){
    double x,y;
    gr->GetPoint(i,x,y);
    cout<<i<<" "<<x<<" "<<y<<endl;
    gr->SetPoint(i,x,y*scale);
  }
}

int LowPassRampPlot(){
  TCanvas *c = new TCanvas("c","c",0,0,1600,900);
  c->Divide(3,2);
  c->cd(1);
  TGraph *gr1 = new TGraph("~/Downloads/FeFoilNoAn2Aux1LowPassFilter031621.txt","%*lg, %lg, %lg");
  scale(gr1);
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerSize(0.25);
  gr1->SetTitle("SBU Fe Foil 99.99%? Purity 1F");
  TGraph *gr2 = new TGraph("~/Downloads/FeFoilNoAn2Aux1RampUp2ndTimeTo4.5TLowPassFilter031621.txt","%*lg, %lg, %lg");
  scale(gr2);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetMarkerSize(0.25);
  gr1->SetTitle("1F Signal (SBU Fe Foil 99.99%? Purity)");
  gr1->Draw("alp");
  gr2->Draw("samelp");
  //gr1->GetXaxis()->SetRangeUser(1.8,4.8);
  //gr1->GetYaxis()->SetRangeUser(0.8,1.03);
  gr1->GetXaxis()->SetTitle("B-field (T)");
  gr1->GetYaxis()->SetTitle("1F (mV)");
  TLegend *lg = new TLegend(0.65,0.2,0.89,0.3);
  lg->AddEntry(gr1,"1st Ramp Up","lp");
  lg->AddEntry(gr2,"2nd Ramp Up","lp");
  lg->Draw();
  c->cd(2);
  TGraph *gr3 = new TGraph("~/Downloads/FeFoilNoAn2Aux1LowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %lg");
  //norm(gr3);
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerColor(kRed);
  gr3->SetLineColor(kRed);
  gr3->SetMarkerSize(0.25);
  gr3->SetTitle("Aux1 DC Signal (SBU Fe Foil 99.99%? Purity)");
  TGraph *gr4 = new TGraph("~/Downloads/FeFoilNoAn2Aux1RampUp2ndTimeTo4.5TLowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %lg");
  //norm(gr4);
  gr4->SetMarkerStyle(8);
  gr4->SetMarkerColor(kBlue);
  gr4->SetLineColor(kBlue);
  gr4->SetMarkerSize(0.25);
  gr3->Draw("alp");
  gr4->Draw("samelp");
  //gr3->GetXaxis()->SetRangeUser(1.8,4.8);
  gr3->GetYaxis()->SetRangeUser(3.9,4.2);
  gr3->GetXaxis()->SetTitle("B-field (T)");
  gr3->GetYaxis()->SetTitle("Aux1 DC (V)");
  lg->Draw();
  c->cd(4);
  TGraph *gr5 = new TGraph("~/Downloads/FeFoilNoAn2Aux1LowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %lg");
  scale(gr5);
  gr5->SetMarkerStyle(8);
  gr5->SetMarkerColor(kRed);
  gr5->SetLineColor(kRed);
  gr5->SetMarkerSize(0.25);
  gr5->SetTitle("2F Signal (SBU Fe Foil 99.99%? Purity)");
  TGraph *gr6 = new TGraph("~/Downloads/FeFoilNoAn2Aux1RampUp2ndTimeTo4.5TLowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %lg");
  scale(gr6);
  gr6->SetMarkerStyle(8);
  gr6->SetMarkerColor(kBlue);
  gr6->SetLineColor(kBlue);
  gr6->SetMarkerSize(0.25);
  gr5->Draw("alp");
  gr6->Draw("samelp");
  //gr3->GetXaxis()->SetRangeUser(1.8,4.8);
  gr5->GetYaxis()->SetRangeUser(1.05,1.45);
  gr5->GetXaxis()->SetTitle("B-field (T)");
  gr5->GetYaxis()->SetTitle("2F (mV)");
  lg->Draw();

  c->cd(5);
  TGraph *gr7 = new TGraph("~/Downloads/FeFoilNoAn2Aux1LowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  //scale(gr7);
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerColor(kRed);
  gr7->SetLineColor(kRed);
  gr7->SetMarkerSize(0.25);
  gr7->SetTitle("Ratio 1F/2F Signal (SBU Fe Foil 99.99%? Purity)");
  TGraph *gr8 = new TGraph("~/Downloads/FeFoilNoAn2Aux1RampUp2ndTimeTo4.5TLowPassFilter031621.txt","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  //scale(gr8);
  gr8->SetMarkerStyle(8);
  gr8->SetMarkerColor(kBlue);
  gr8->SetLineColor(kBlue);
  gr8->SetMarkerSize(0.25);
  gr7->Draw("alp");
  gr8->Draw("samelp");
  gr7->GetXaxis()->SetTitle("B-field (T)");
  gr7->GetYaxis()->SetTitle("Ratio 1F/2F");
  lg->Draw();
 c->cd(3);
 TGraph *gr9 = (TGraph*)gr1->Clone("gr9");
 TGraph *gr10 = (TGraph*)gr2->Clone("gr10");
 for(int i=0;i<gr9->GetN();++i){
   double x,y,x1,y1;
   gr9->GetPoint(i,x,y);
   gr3->GetPoint(i,x1,y1);
   gr9->SetPoint(i,x,y/y1);
 }
 for(int i=0;i<gr10->GetN();++i){
   double x,y,x1,y1;
   gr10->GetPoint(i,x,y);
   gr4->GetPoint(i,x1,y1);
   gr10->SetPoint(i,x,y/y1);
 }
 gr9->Draw("alp");
 gr9->SetTitle("Ratio 1F/Aux1 DC Signal");
 gr9->GetYaxis()->SetTitle("Ratio 1F/Aux1 DC");
 gr10->Draw("samelp");
 lg->Draw();
 c->cd(6);
 TGraph *gr11 = (TGraph*)gr9->Clone("gr11");
 gr11->SetTitle("Ratio 1F/Aux1 DC Signal Zoomed");
 TGraph *gr12 = (TGraph*)gr10->Clone("gr12");
 gr11->GetXaxis()->SetRangeUser(1.8,4.8);
 gr11->GetYaxis()->SetRangeUser(0.68,0.77);
 gr11->Draw("alp");
 gr12->Draw("samelp");
 c->SaveAs("LowPassRampUp.png");
  return 0;
}
