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
int MagFoilAngle(){
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  TF1 *f = new TF1("f",saturation,0,4,2);
  // f->Draw();
  f->SetParameters(89,0.052);
  f->SetLineColor(kRed);
  TF1 *f1 = new TF1("f1",saturation,0,4,2);
  f1->SetParameters(89,0.052);
  double x[5] = {1,2.4,2.8,3.2,4}, xe[5] = {0,0,0,0,0};
  double y[5] = {0.0246,0.05207,0.05269,0.05227,0.05243}, ye[5] = {0.00015,0.00016,0.00011,0.0001,0.0001};
  double azz_sim[5] = {0.752561,0.752765,0.752909,0.753319,0.753826};
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  TGraphErrors *gr = new TGraphErrors(5,x,y,xe,ye);
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kBlue);
  gr->SetMarkerStyle(8);
  gr->SetTitle("Moller Asymmetry versus Field");
  gStyle->SetOptFit(1111);
  gr->Fit(f);
  gr->Draw("ap");
  gr->GetYaxis()->SetTitle("A_{meas}");
  gr->GetXaxis()->SetTitle("Target Coil Field (T)");
  gPad->Update();
  c->SaveAs("FoilSaturationCurve.png");
  TCanvas *c1 = new TCanvas("c1","c1",800,0,1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  TGraphErrors *grc = new TGraphErrors();
  TGraphErrors *grcu = new TGraphErrors();
  int n=0;
  for(int i=0;i<5;++i){
    grc->SetPoint(i,x[i],y[i]*azz_sim[4]/azz_sim[i]);
    grc->SetPointError(i,0,ye[i]*azz_sim[4]/azz_sim[i]);
    if(i>2 ){
      grcu->SetPoint(n,x[i],y[i]*azz_sim[4]/azz_sim[i]);
      grcu->SetPointError(n,0,ye[i]*azz_sim[4]/azz_sim[i]);
      ++n;
    }
  }
  grc->SetMarkerColor(kBlue);
  grc->SetLineColor(kBlue);
  grc->SetMarkerStyle(8);
  grc->SetTitle("Moller Asymmetry Corrected for A_{zz} Evolution vs. Field");
  gStyle->SetOptFit(1111);
  f->SetRange(2.3,2.9);
  grc->Fit(f,"r");
  grc->Draw("ap");
  grc->GetYaxis()->SetTitle("A_{meas}");
  grc->GetXaxis()->SetTitle("Target Coil Field (T)");
  c1->cd(2);
  grcu->SetMarkerColor(kRed);
  grcu->SetLineColor(kRed);
  grcu->SetMarkerStyle(8);
  f1->SetLineColor(kRed);
  f1->SetParameters(89.8,0.53);
  f1->SetRange(2,5);
  grcu->Fit(f1,"r");
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Moller Asymmetry Corrected for A_{zz} Evolution vs. Field");
  mg->Add(grc);
  mg->Add(grcu);
  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("A_{meas}");
  mg->GetXaxis()->SetRangeUser(2,4.5);
  mg->GetYaxis()->SetRangeUser(0.051,0.054);
  mg->GetXaxis()->SetTitle("Target Coil Field (T)");
  f->SetRange(2,4.2);
  f->Draw("same");
  cout<<"Blue curve of lower 2 pts. evaluated at 4.0 T: "<<f->Eval(4)<<endl;
  gPad->Update();
  c1->SaveAs("FoilSaturationCurve_AzzCorrected.png");
  double low=87.76, high=88.94;
  f->SetParameters(low,1.0);
  f1->SetParameters(high,1.0);
  cout<<f->Eval(4.0)<<endl;
  TPaveText *pt = new TPaveText(0.5,0.6,0.93,0.7,"ndc");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("Top 2 pts scaled from %0.2f#circ to %0.2f#circ",low, high));
  //pt->AddText(Form("Extra 0.3%% point-to-point error added in quadrature"));
  if(1){
    cout<<"Scaling graph\n";
    double scale[5] = {1,1,1,f1->Eval(3.5)/f->Eval(3.5),f1->Eval(4)/f->Eval(4)};
    //    cout<<scale[3]<<" "<<scale[4]<<endl;
    f->SetParameters(89,0.052);
    f->SetRange(2,5);
    for(int i=0;i<5;++i){
      grc->SetPoint(i,x[i],scale[i]*y[i]*azz_sim[4]/azz_sim[i]);
      grc->SetPointError(i,0,scale[i]*ye[i]*azz_sim[4]/azz_sim[i]);
      //grc->SetPointError(i,0,sqrt(pow(scale[i]*ye[i]*azz_sim[4]/azz_sim[i],2)+pow(0.000158,2)));
    }
    grc->Fit(f,"r");
    grc->Draw("ap");
    grc->GetXaxis()->SetRangeUser(2,4.5);
    grc->GetYaxis()->SetRangeUser(0.051,0.054);
    pt->Draw();
    gPad->Update();
  }
  // TCanvas *cp = new TCanvas("cp","cp",0,0,800,600);
  // double yp[2] = {(0.054319+0.055147)/2.0,(0.054515+0.055180)/2.0};
  // double ype[2] = {0.00017/sqrt(2), 0.00017/sqrt(2)};
  // double xp[2] = {3.5,4};
  // TGraphErrors *grp = new TGraphErrors(2,xp,yp,xe,ype);
  // grp->SetMarkerColor(kRed);
  // grp->SetLineColor(kRed);
  // grp->SetMarkerStyle(8);
  // grp->Draw("ap");
  // f->SetParameters(88,0.54);
  // grp->Fit(f);
  TCanvas *cp = new TCanvas("cp","cp",0,0,800,600);
  TGraph *gr1 = new TGraph("kerr.dat","%lg %*lg %lg");
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerStyle(8);
  gr1->Draw("ap");
  f->SetParameters(85,1);
  gr1->Fit(f); 

  TGraph *gr2 = new TGraph("kerr.dat","%lg %*lg%*lg %lg");
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetMarkerStyle(8);
  gr2->Draw("samep");
  f->SetLineColor(kBlue);
  f->SetParameters(85,1);
  gr2->Fit(f);
  gr1->SetTitle("Kerr Data LockIn Signal (SBU=Blue) (JLab=Red)");
  gr1->GetXaxis()->SetTitle("Applied Field (B)");
  gr1->GetYaxis()->SetTitle("LockIn Signal Ratio: F1/F2");
  //gr1->GetXaxis()->SetRangeUser(2,4.5);
 return 0;
}
