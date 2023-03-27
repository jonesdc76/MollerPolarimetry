#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TRandom.h>
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
  return s;
}

int FoilNormalScan(){
  TCanvas *c = new TCanvas("c","c",0,0,1400,600);
  c->Divide(2,1);
  c->cd(1);
  double x[5] = {86,88,90,92,94}, xe[5]={0,0,0,0,0}, xalt[5]={-4,-2,0,2,4};
  double y[5] = {4.90421,5.07748,5.10825,5.04485,4.92099}, ye[5]={0.01459,0.03781,0.0085,0.01071,0.01152};
  TGraphErrors *gre = new TGraphErrors(5,xalt,y,xe,ye);
  gre->SetMarkerStyle(8);
  //gre->Draw("ap");
  double par[2] = {0,y[2]};
  TGraph *gr = new TGraph();
  double xt = 2.58;
  for(int i=0;i<90;++i){
    double xtmp = double(85.5+0.1*i);
    par[0] = xtmp <= 90 ? xtmp : 180-xtmp;
    par[1] = 5.10825;
    double ytmp = 5.10825*saturation(&xt, par);
    xtmp -= 90;
    gr->SetPoint(i,xtmp,ytmp);
    printf("%f   %f\n",xtmp,ytmp);
  }
  gStyle->SetStatW(0.15);
  gr->SetLineColor(kGreen+2);
  gr->SetLineWidth(3);
  gr->Draw("al");
  gr->SetTitle("Moller Asymmetry vs Foil Angle");
  gr->GetXaxis()->SetTitle("Moller Asymmetry (%)");
  gr->GetXaxis()->SetTitle("Foil Angle (deg)");
  gr->GetYaxis()->SetLimits(4.85,5.2);
  gr->GetYaxis()->SetRangeUser(4.85,5.2);
  gPad->Update();
  gre->Draw("samep");
  gPad->Update();
  TLegend *tl = new TLegend(0.1,0.93,0.43,0.75);
  tl->AddEntry(gre,"Data scan","lp");
  tl->AddEntry(gr,"Stoner-Wohlfarth Model","l");
  TF1 *f = new TF1("f","pol2",0,1);
  gStyle->SetOptFit(1111);
  TFitResultPtr frpt = gre->Fit(f,"s");
  f->SetLineWidth(3);
  f->SetLineStyle(7);
  tl->AddEntry(f,"Quadratic Fit","l");
  tl->Draw();
  TPaveText *pt = new TPaveText(0.255,0.15,0.75,0.25,"ndc");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  double b = f->GetParameter(1), a = f->GetParameter(2);
  double be = f->GetParError(1), ae = f->GetParError(2);
  double var = pow(be/(2*a),2);
  var += pow(b*ae/(2*a*a),2);
  cout<<"Var "<<var<< " "<<2*pow(frpt->GetCovarianceMatrix()[1][2],2)*(-1/(2*a))*(b/(2*a*a))<<endl;
  var += 2*pow(frpt->GetCovarianceMatrix()[1][2],2)*(-1/(2*a))*(b/(2*a*a));
  pt->AddText(Form("Maximum at %0.2f#circ+/-%0.2f#circ\n",-b/(2*a),sqrt(var)));
  pt->SetTextColor(kRed);
  pt->Draw();
  cout<<frpt->GetCovarianceMatrix()[0][0]<<endl;
  cout<<frpt->GetCovarianceMatrix()[0][1]<<endl;
  cout<<frpt->GetCovarianceMatrix()[1][2]<<endl;
  cout<<frpt->GetCovarianceMatrix()[1][1]<<endl;
  cout<<frpt->GetCovarianceMatrix()[2][2]<<endl;
  TGraphErrors *grt;
  TRandom r(0);
  TH1D *h = new TH1D("h","h",100,89,91);
  double yt[5] = {0,0,0,0,0};
  for(int i=0;i<1000;++i){
    for(int j=0;j<5;++j)
      yt[j] = r.Gaus(y[j],ye[j]);
    grt = new TGraphErrors(5,x,yt,xe,ye);
    grt->Fit(f,"q");
    h->Fill(-f->GetParameter(1)/(2*f->GetParameter(2)));
  }
  c->cd(2);
  h->Draw();
  xt = 3.2;
  par[0] = 89.75;
  cout<<saturation(&xt,par)*100<<"%"<<endl;
  return 0;
}
