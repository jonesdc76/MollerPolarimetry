#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
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
int MagFoilAngle1(){
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  TF1 *f = new TF1("f",saturation,1.8,4.5,2);
  // f->Draw();
  f->SetParameters(89,2);
  f->SetLineColor(kRed);
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  TGraph *gr = new TGraph("FeFoil4.5to0TFoilAngleNeg2.535NoA2_031021.txt","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  gr->SetTitle("1F/2F Signal versus Field");
  if(0){
    //subtract 0.1% per Tesla for high field magnetization not included
    //in saturation model
    double x, y;
    for(int i=0;i<gr->GetN();++i){
      gr->GetPoint(i,x,y);
      gr->SetPoint(i,x,y*(1-0.001*x));
    }
    gr->SetTitle("1F/2F Signal versus Field (Corrected 0.1%/T for High Field Magnetization)");

  }
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kBlue);
  gr->SetMarkerStyle(8);
  gStyle->SetOptFit(1111);
  gr->Fit(f,"r");
  gr->Draw("ap");
  gr->GetYaxis()->SetTitle("1F/2F");
  gr->GetXaxis()->SetTitle("Target Coil Field (T)");

 return 0;
}
