#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>
#include <TCanvas.h>
const double pi = 3.1415927, eps = 0.00001;

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
  //get within 0.1%
  double s = 0, s1=0.25, s2 = 0.75;
  for(int i=0;i<10;++i){
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
  //cout<<f1<<" "<<f2<<" "<<s1<<" "<<s2<<" "<<s<<endl;
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
int MagFoilAngle(){
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  TF1 *f = new TF1("f",saturation,0,4,2);
  // f->Draw();
  f->SetParameters(89,0.052);
  double x[5] = {1,2.4,2.8,3.2,4}, xe[5] = {0,0,0,0,0};
  double y[5] = {0.0246,0.05207,0.05269,0.05227,0.05243}, ye[5] = {0.00015,0.00016,0.00011,0.0001,0.0001};
  double azz_sim[5] = {0.752561,0.752765,0.752909,0.753319,0.753826};
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  TGraphErrors *gr = new TGraphErrors(5,x,y,xe,ye);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(8);
  gr->SetTitle("Moller Asymmetry versus Field");
  gStyle->SetOptFit(1111);
  gr->Fit(f);
  gr->Draw("ap");
  gr->GetYaxis()->SetTitle("A_{meas}");
  gr->GetXaxis()->SetTitle("Target Coil Field (T)");
  gPad->Update();
  c->SaveAs("FoilSaturationCurve.png");
  TCanvas *c1 = new TCanvas("c1","c1",800,0,800,600);
  TGraphErrors *grc = new TGraphErrors();
  double scale[5] = {1,1,1,1.00,1.00};
  for(int i=0;i<5;++i){
    grc->SetPoint(i,x[i],scale[i]*y[i]*azz_sim[4]/azz_sim[i]);
    grc->SetPointError(i,0,scale[i]*ye[i]*azz_sim[4]/azz_sim[i]);
  }
  grc->SetMarkerColor(kBlue);
  grc->SetMarkerStyle(8);
  grc->SetTitle("Moller Asymmetry Corrected for A_{zz} Evolution vs. Field");
  gStyle->SetOptFit(1111);
  grc->Fit(f);
  grc->Draw("ap");
  grc->GetYaxis()->SetTitle("A_{meas}");
  grc->GetXaxis()->SetTitle("Target Coil Field (T)");
  gPad->Update();
  f->SetParameters(89,0.052);
  double x32=f->Eval(3.2), x40=f->Eval(4.0);
  f->SetParameters(87,0.052);
  cout<<f->Eval(3.2)/x32<<"    "<<f->Eval(4.0)/x40<<endl;
  c1->SaveAs("FoilSaturationCurve_AzzCorrected.png");
  return 0;
}
