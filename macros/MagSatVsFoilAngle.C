#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
const double pi = 3.1415927, eps = 0.000001;

//The functional form gives the field as a function of angle and saturation,
//but we need saturation as a function of field and angle.
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
int MagSatVsFoilAngle(){
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  //c->SetGrid();
  const int N=5, N1=1000;
  TString name[N] = {"83 deg", "85 deg", "87 deg","89 deg", "90 deg"};
  int col[N]={kBlack, kGreen+2, kBlue, kRed, kBlack};
  int style[N]={7, 8, 9, 10, 1};
  double x[N][N1], y[N][N1], par[N]={86,87,88,89,90};
  for(int i=0;i<N;++i)name[i] = Form("%0.1f deg",par[i]);
  TGraph *gr[N];
  TLegend *tl = new TLegend(0.726,0.25,0.94,0.5);
  tl->SetBorderSize(0);
  for(int i=0;i<N;++i){
    double p[2]={par[i],1};
    for(int j=0;j<N1;++j){
      x[i][j] = j*0.0042;
      y[i][j] = saturation(&x[i][j], p);
      cout<<x[i][j]<<" "<<y[i][j]<<endl;
    }
    gr[i] = new TGraph(N1, x[i], y[i]);
    gr[i]->SetLineWidth(4);
    gr[i]->SetLineColor(col[i]);
    gr[i]->SetLineStyle(style[i]);
    gr[i]->SetMarkerColor(col[i]);
    //gr[i]->SetMarkerStyle(8);
    if(i==0){
      gr[i]->SetTitle("");//Degree of Magnetic Saturation versus Foil Angle");
      gr[i]->Draw("al");
      gr[i]->GetXaxis()->SetTitle("Applied Field (T)");
      gr[i]->GetXaxis()->SetTitleSize(0.044);
      gr[i]->GetYaxis()->SetTitleSize(0.044);
      gr[i]->GetYaxis()->SetTitle("Fractional Saturation");
      gr[i]->GetXaxis()->SetRangeUser(2,4.1);
      gr[i]->GetYaxis()->SetRangeUser(0.95,1.004);
    }else{
      gr[i]->Draw("samec");
    }
  }
  for(int i=N;i>0;i--) tl->AddEntry(gr[i-1],name[i-1].Data(),"lp");

  tl->Draw();
  c->SaveAs("StonerCurves.pdf");
 return 0;
}
