#include <iostream>
#include <cstdio>
#include "TGaxis.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMultiGraph.h"

//Returns theoretical magnetization value for given Hi in Oersted, T in Kelvin
//Curves from Pauthenet Mar 1982 "Spin Waves in nickel, iron and yttrium-iron
//garnet" equation 10 and Table 1
//Adjustable parameter is spontaneous magnetization (nominally 58.858 emu) 
double magnetizationNi( double Hi, double T, double spont_m = 58.858,  
		      int nTerms = 1000){
  //calculate 307e-6*T^(3/2)*F(3/2, 1.378*Hi/T) from equation (10) at T and Hi
  const double a3_2 = -154e-6, a5_2 = -70.3e-8, a7_2 = -3.8e-10, b = 1.478e-4;
  double s = 3.0/2.0, term1 = 0, term2 = 0, term3 = 0, term4 = 0;

  for(int i=1;i<nTerms;++i){
    term1 += a3_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  s = 5.0/2.0;
  for(int i=1;i<nTerms;++i){
    term2 += a5_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  s = 7.0/2.0;
  for(int i=1;i<nTerms;++i){
    term3 += a7_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  //my linear parameterization of Chi term using fit to values in Table 1
  double intercept = 1.60986e-6, slope = -1.83031e-9;
  term4 = (intercept + slope * T) * Hi;
  double mag = spont_m + term1 + term2 + term3 + term4;
  return mag;
}

//Returns theoretical magnetization value for given Hi in Oersted, T in Kelvin
//Curves from Pauthenet Mar 1982 "Spin Waves in nickel, iron and yttrium-iron
//garnet" equation 9 and Table 1
//Adjustable parameter is spontaneous magnetization (nominally 222.678 emu) 
double magnetizationFe( double Hi, double T, double spont_m = 222.678,
		      int nTerms = 10000){
  gStyle->SetTitleW(0.8);
  //calculate 307e-6*T^(3/2)*F(3/2, 1.378*Hi/T) from equation (9) at T and Hi
  const double a3_2 = -307e-6, a5_2 = -22.8e-8, b = 1.378e-4;
  double s = 3.0/2.0, term1 = 0, term2 = 0, term3 = 0;
  for(int i=1;i<nTerms;++i){
    term1 += a3_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  s = 5.0/2.0;
  for(int i=1;i<nTerms;++i){
    term2 += a5_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  //my linear parameterization of Chi term using fit to values in Table 1
  double intercept = 3.644e-6, slope = 5.0434e-10;
  term3 = (intercept + slope * T) * Hi;
  double mag = spont_m + term1 + term2 + term3;
  return mag;
}
const int N=30;
void HeatingCorrections(double startT = 294, bool fractional = 1){
  double H_app_Ni = 20000, H_app_Fe = 40000;//2T Ni and 4 T Fe applied field 
  const double H_sat_Ni = 6179, H_sat_Fe = 21570;
  TMultiGraph *mg = new TMultiGraph();
  double x[N], yNi[N], yFe[N];
  double Mroom_Ni = magnetizationNi(H_app_Ni-H_sat_Ni, startT);
  double Mroom_Fe = magnetizationFe(H_app_Fe-H_sat_Fe, startT);
  for(int i=0;i<N;++i){
    x[i] = i+startT;
    yNi[i] = magnetizationNi(H_app_Ni-H_sat_Ni, x[i]);
    yFe[i] = magnetizationFe(H_app_Fe-H_sat_Fe, x[i]);
    if(fractional){
      yNi[i]/=Mroom_Ni;
      yFe[i]/=Mroom_Fe;
    }else{
      yNi[i]-=Mroom_Ni;
      yFe[i]-=Mroom_Fe;
    }
  }
  TGraph *grNi = new TGraph(N,x,yNi);
  grNi->SetLineColor(kRed);
  grNi->SetLineWidth(2);
  grNi->SetMarkerColor(kRed);
  grNi->SetMarkerStyle(4);
  mg->Add(grNi);
  TGraph *grFe = new TGraph(N,x,yFe);
  grFe->SetLineColor(kBlack);
  grFe->SetLineWidth(2);
  grFe->SetMarkerColor(kBlack);
  grFe->SetMarkerStyle(4);
  mg->Add(grFe);
  TString title = "Magnetization Correction vs Temperature for Fe and Ni Foils";
  if(fractional) title = "Reduced Magnetization vs Temperature for Fe and Ni Foils";
  mg->SetTitle(title.Data());
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("T (K)");
  mg->GetYaxis()->SetTitle("Magnetization Correction (emu/g)");
  if(fractional){
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    mg->GetYaxis()->SetTitle("#frac{M(T)}{M(294 K)}");
    //mg->Draw("ap");
  }
  // mg->GetXaxis()->SetTitleOffset(1.2);
  // mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->SetLimits(-0.8, 0.1);
  mg->GetYaxis()->SetRangeUser(-0.8, 0.1);
  if(fractional){
    mg->GetYaxis()->SetLimits(0.985, 1.003);
    mg->GetYaxis()->SetRangeUser(0.985, 1.003);
  }
  mg->Draw("ap");
  TF1 *fFe = new TF1("fFe","pol1",294.3,323.5);
  fFe->SetLineColor(kBlack);
  fFe->SetLineWidth(3);
  grFe->Fit(fFe,"R");
  TF1 *fNi = new TF1("fNi","pol1",294.3,323.5);
  fNi->SetLineColor(kRed);
  fNi->SetLineWidth(3);
  grNi->Fit(fNi,"R");
  TLegend *tl = new TLegend(0.6,0.7,0.89,0.89);
  tl->SetBorderSize(0);
  if(fractional){
    tl->AddEntry(grFe,Form("Fe: %0.3f (%%/K)", fFe->GetParameter(1)*100),"lp");    
    tl->AddEntry(grNi,Form("Ni: %0.3f (%%/K)",fNi->GetParameter(1)*100),"lp");
  }else{
    tl->AddEntry(grFe,Form("Fe: %0.3f (emu/g#circC)", fFe->GetParameter(1)),"lp");
    tl->AddEntry(grNi,Form("Ni: %0.3f (emu/g#circC)",fNi->GetParameter(1)),"lp");
  }
  tl->Draw();
  c->SaveAs("target_heating_correction.pdf");
  c->SaveAs("../nim/figures/target_heating_correction.pdf");
  return;
}
