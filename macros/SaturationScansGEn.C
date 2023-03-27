#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TCanvas.h>
const double pi = 3.1415927, eps = 0.000001;

//The functional form gives the field as a function of angle and saturation,
//but we need satuation as a function of field and angle.
//This function solves the equation to give saturation.
//Arg 1, x, is a pointer to the field strength variable in Tesla
//Arg 2, par, is an array of 2 fit parameters:
//       par[0] = foil angle in degrees
//       par[1] = normalization.
double saturation(double *x, double *par){
  double f = x[0];
  double a = par[0]*pi/180.0;
  if(par[0] == 90){
    return (f>2.2 ? 1.0 : f/2.2);
  }
  if(f<0) return 0;
  //start by getting within 1 ppm
  double s = 0, s1=0.25, s2 = 0.75;
  for(int i=0;i<20;++i){
    if(abs(f+sin(2*(acos(s1)-a))/sin(acos(s1))*1.1) < 
       abs(f+sin(2*(acos(s2)-a))/sin(acos(s2))*1.1)){
      s2 = s1 + pow(0.5,i+3); 
      s1 = s1 - pow(0.5,i+3); 
    }else{
      s1 = s2 - pow(0.5,i+3); 
      s2 = s2 + pow(0.5,i+3); 
    }
  }
  double f1 = -sin(2*(acos(s1)-a))/sin(acos(s1))*1.1;
  double f2 = -sin(2*(acos(s2)-a))/sin(acos(s2))*1.1;
  s = (f-f1)*(s2-s1)/(f2-f1)+s1;
  //printf("%0.12f %0.12f %0.12f %0.12f %0.12f\n",f1,f2,s1,s2,s);
  while(abs(f+sin(2*(acos(s)-a))/sin(acos(s))*1.1)>eps){
    s1 = s2; s2 = s;
    //cout<<s1<<" "<<s2;
    f1 = -sin(2*(acos(s1)-a))/sin(acos(s1))*1.1;
    f2 = -sin(2*(acos(s2)-a))/sin(acos(s2))*1.1;
    s = (f-f1)*(s2-s1)/(f2-f1)+s1;
    //cout<<" "<<s<<" "<<f1<<" "<<f2<<endl;
  }
  return s*par[1];
}

//Uses saturation function to fit Moller asymmetry vs field data
//to extract the best fit foil angle and normalization
int SaturationScansGEn(){
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  gStyle->SetFitFormat("6.4g");
  TPaveText *ptt2 = new TPaveText(0.6,0.45,0.94,0.65,"ndc");
  TPaveText *ptt3 = new TPaveText(0.6,0.2,0.94,0.4,"ndc");
  ptt2->SetLineStyle(10);
  ptt2->SetLineColor(kGreen+2);
  ptt2->SetShadowColor(0);
  ptt2->SetLineWidth(3);
  ptt2->SetFillColor(0);
  ptt3->SetLineStyle(2);
  ptt3->SetLineWidth(3);
  ptt3->SetShadowColor(0);
  ptt3->SetFillColor(0);
  TF1 *f = new TF1("f",saturation,0,4,2);
  // f->Draw();
  f->SetParNames("Foil Angle", "Normalization");
  f->SetParameters(89,0.052);
  TCanvas *c = new TCanvas("c","c",0,0,1000,700);
  TMultiGraph *mg = new TMultiGraph();
  double x[5]={4.02,2.8,2.5,2.25,2.0},xe[5]={0,0,0,0,0};
  double y[5]={82.83,82.74,82.86,81.64,76.42},ye[5]={0.22,0.21,0.21,0.2,0.29};
  TGraphErrors *gr2 = new TGraphErrors(5,x,y,xe,ye);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.4);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(kGreen+2);
  gr2->SetLineColor(kGreen+2);
  //gr2->Draw("ap");

  double x3[7]={4.02,3.2,2.8,2.5,2.25,2.35,2.0},x3e[7]={0,0,0,0,0,0,0};
  double y3[7]={82.47,82.76,82.56,82.05,79.76,80.67,74.64},y3e[7]={0.15,0.21,0.2,0.2,0.2,0.25,0.31};
  TGraphErrors *gr3 = new TGraphErrors(7,x3,y3,x3e,y3e);
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerSize(1.4);
  gr3->SetLineWidth(2);
  gr3->SetMarkerColor(kBlack);
  gr3->SetLineColor(kBlack);
  //gr3->Draw("samep");
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Draw("alp");
  mg->SetTitle("Polarization versus Target Field");
  mg->GetXaxis()->SetTitle("Target Field (T)");
  mg->GetYaxis()->SetTitle("Polarization(%)");
  mg->Draw("ap");

  TLegend *lg = new TLegend(0.28,0.2,0.48,0.35);
  lg->SetBorderSize(0);
  lg->AddEntry(gr2,"Foil 2","lp");
  lg->AddEntry(gr3,"Foil 3","lp");
  gStyle->SetOptFit(0);
  f->SetRange(2.0,4.1);
  f->SetLineWidth(3);
  f->SetLineColor(kGreen+2);
  f->SetLineStyle(10);
  f->SetParameters(89,83);
  gr2->Fit(f,"r");
  ptt2->SetTextColor(kGreen+2);
  ptt2->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  ptt2->AddText(Form("P-value                %0.4f",f->GetProb()));
  ptt2->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f",f->GetParameter(0), f->GetParError(0)));
  ptt2->AddText(Form("Max Polarization (%%)   %0.2f #pm %0.2f",f->GetParameter(1), f->GetParError(1)));
  f->SetLineStyle(8);
  f->SetLineColor(kBlack);
  gr3->Fit(f,"r");
  ptt3->SetTextColor(kBlack);
  ptt3->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  ptt3->AddText(Form("P-value                %0.4f",f->GetProb()));
  ptt3->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f",f->GetParameter(0), f->GetParError(0)));
  ptt3->AddText(Form("Max Polarization (%%)   %0.2f #pm %0.2f",f->GetParameter(1), f->GetParError(1)));
  mg->Draw("ap");
  lg->Draw();
  ptt2->SetTextAlign(11);
  ptt3->SetTextAlign(11);
  ptt3->Draw("l");
  gPad->Update();
  ptt2->Draw();
  ptt3->Draw();
  //  c->SaveAs("FoilSaturationCurve.pdf");
  return 0;

}
