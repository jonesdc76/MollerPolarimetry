#include <iostream>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLine.h>
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
int SaturationScansVsAngleGEn(){
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.35);
  gStyle->SetFitFormat("6.4g");

  int col[5] = {kGreen+2,kRed,kBlue,kBlack,kMagenta};
  int style[5] = {21,22,23,24,25};
  TF1 *f = new TF1("f",saturation,0,4,2);
  // f->Draw();
  f->SetParNames("Foil Angle", "Normalization");
  f->SetParameters(89,5.2);
  gStyle->SetOptFit(0);
  f->SetRange(2.0,4.1);
  f->SetLineWidth(3);
  f->SetLineStyle(10);
  f->SetParameters(89,5.1297175);
  //f->FixParameter(1,5.1297175);
  TCanvas *c = new TCanvas("c","c",0,0,1000,700);
  TMultiGraph *mg = new TMultiGraph();
  double x[5]={2.3,2.4,2.5,3.2,4.0},xe[5]={0,0,0,0,0};
  //for(int i=0;i<5;++i)x[i]+=0.08;
  double y0[5]={5.0809,5.0970,5.1177,5.13415,5.1157}, y0e[5]={0.0072,0.0052,0.0052,0.00764,0.0048};
  double y1[5]={5.0222,5.0549,5.091,5.1395,5.1168},y1e[5]={0.0068,0.0066,0.0067,0.0050,0.0048};
  double y2[5]={4.9400,4.9835,5.0592,5.1385,5.1231},y2e[5]={0.0235,0.0074,0.0072,0.0052,0.0052};
  double y4[5]={4.8142,4.8958,4.9462,5.1100,5.1266},y4e[5]={0.0099,0.0077,0.0069,0.0052,0.0053};

  int n=0;
  TPaveText *pt = new TPaveText(0.7,0.15,0.94,0.75,"ndc");
  pt->SetShadowColor(0);
  //  pt->SetLineWidth(0);
  pt->SetFillColor(0);
  pt->SetTextAlign(11); 
  TGraphErrors *gr0 = new TGraphErrors(5,x,y0,xe,y0e);
  f->SetLineColor(col[n]);
  f->SetLineStyle(8);
  f->SetLineWidth(3);
  gr0->Fit(f);

  TText *tt;TLine *tl;
  pt->SetTextSize(0.02);
  //pt->AddText(" ");
  tt = pt->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("P-value                %0.4f\n",f->GetProb()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f\n",f->GetParameter(0),f->GetParError(0)));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Max Asymmetry (%%)   %0.4f #pm %0.4f",f->GetParameter(1),f->GetParError(1)));
  tt->SetTextColor(col[n]);
  pt->SetTextSize(0.0005);
  pt->AddText(" ");
  pt->SetTextSize(0.02);
  tl = pt->AddLine();
  tl->SetLineColor(col[n]); 
  gr0->SetMarkerStyle(style[n]);
  gr0->SetMarkerSize(1.4);
  gr0->SetLineWidth(2);
  gr0->SetMarkerColor(col[n]);
  gr0->SetLineColor(col[n]);

  ++n;
  TGraphErrors *gr1 = new TGraphErrors(5,x,y1,xe,y1e);
  f->SetLineColor(col[n]);
  f->SetLineColor(col[n]);
  f->SetLineStyle(8);
  f->SetLineWidth(3);
  gr1->Fit(f);
  tt = pt->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("P-value                %0.4f\n",f->GetProb()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f\n",f->GetParameter(0),f->GetParError(0)));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Max Asymmetry (%%)   %0.4f #pm %0.4f",f->GetParameter(1),f->GetParError(1)));
  tt->SetTextColor(col[n]);
  pt->SetTextSize(0.005);
  pt->AddText(" ");
  pt->SetTextSize(0.02);
  tl = pt->AddLine();
  tl->SetLineColor(col[n]); 

  gr1->Fit(f);
  gr1->SetMarkerStyle(style[n]);
  gr1->SetMarkerSize(1.4);
  gr1->SetLineWidth(2);
  gr1->SetMarkerColor(col[n]);
  gr1->SetLineColor(col[n]);
  ++n;
  TGraphErrors *gr2 = new TGraphErrors(5,x,y2,xe,y2e);
  f->SetLineColor(col[n]);
  f->SetLineColor(col[n]);
  f->SetLineStyle(8);
  f->SetLineWidth(3);
  gr2->Fit(f);
  tt = pt->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("P-value                %0.4f\n",f->GetProb()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f\n",f->GetParameter(0),f->GetParError(0)));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Max Asymmetry (%%)   %0.4f #pm %0.4f",f->GetParameter(1),f->GetParError(1)));
  tt->SetTextColor(col[n]);
  pt->SetTextSize(0.005);
  pt->AddText(" ");
  pt->SetTextSize(0.02);
  tl = pt->AddLine();
  tl->SetLineColor(col[n]); 

  gr2->SetMarkerStyle(style[n]);
  gr2->SetMarkerSize(1.4);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(col[n]);
  gr2->SetLineColor(col[n]);
  ++n;
  TGraphErrors *gr4 = new TGraphErrors(5,x,y4,xe,y4e);
  f->SetLineColor(col[n]);
  f->SetLineColor(col[n]);
  f->SetLineStyle(8);
  f->SetLineWidth(3);
  gr4->Fit(f);
  tt = pt->AddText(Form("#chi^{2}/NDF                %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("P-value                %0.4f\n",f->GetProb()));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Foil Angle (degrees)   %0.2f #pm %0.2f\n",f->GetParameter(0),f->GetParError(0)));
  tt->SetTextColor(col[n]);
  tt = pt->AddText(Form("Max Asymmetry (%%)   %0.4f #pm %0.4f",f->GetParameter(1),f->GetParError(1)));
  tt->SetTextColor(col[n]);

  gr4->SetMarkerStyle(style[n]);
  gr4->SetMarkerSize(1.4);
  gr4->SetLineWidth(2);
  gr4->SetMarkerColor(col[n]);
  gr4->SetLineColor(col[n]);

  mg->Add(gr0);
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr4);
  mg->Draw("alp");
  mg->SetTitle("Measured Asymmetry versus Target Field");
  mg->GetXaxis()->SetTitle("Target Field (T)");
  mg->GetYaxis()->SetTitle("Asymmetry (%)");
  mg->Draw("ap");

  TLegend *lg = new TLegend(0.28,0.2,0.48,0.35);
  lg->SetBorderSize(0);
  lg->AddEntry(gr0,"0 degrees","lp");
  lg->AddEntry(gr1,"1 degrees","lp");
  lg->AddEntry(gr2,"2 degrees","lp");
  lg->AddEntry(gr4,"4 degrees","lp");
  mg->Draw("ap");
  lg->Draw();
  pt->Draw(); 

  //  c->SaveAs("FoilSaturationCurve.pdf");
  return 0;

}
