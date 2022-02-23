/*Plots the saturation magnetization data for Ni from numerous papers including
Weiss 1929 "La Saturation Absolue Des Ferromagnetiques"

Danan 1959 "On the interpretation of the magnetization measurements of pure polycrystalline iron and nickel in the vicinity of saturation"

Araj et al 1967 "Electrical Resistivity Studies of Chromium-Rich Chromium-Cobalt Alloys

Crangle et al 1971 "The Magnetization of Pure Iron and Nickel"

Shull 2000 "Absolute magnetic moment measurements of nickel spheres"

Plots are converted where possible to be versus applied field H

 */
#include <iostream>
#include <cstdio>
#include "TGaxis.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLine.h"

//Returns theoretical magnetization value for given Hi in kOe, T in Kelvin
//Curves from Pauthenet Mar 1982 "Spin Waves in nickel, iron and yttrium-iron
//garnet" equation 10 and Table 1
//Adjustable parameter is spontaneous magnetization (nominally 58.858 emu) 
double magnetization( double Hi, double T, double spont_m = 58.858,  
		      int nTerms = 1000){
  //calculate 307e-6*T^(3/2)*F(3/2, 1.378*Hi/T) from equation (10) at T and Hi
  const double a3_2 = -154e-6, a5_2 = -70.3e-8, a7_2 = -3.8e-10, b = 1.478e-4;
  double s = 3.0/2.0, term1 = 0, term2 = 0, term3 = 0, term4 = 0;

  for(int i=1;i<nTerms;++i){
    term1 += a3_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi*1000.0/T);
  }
  s = 5.0/2.0;
  for(int i=1;i<nTerms;++i){
    term2 += a5_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi*1000.0/T);
  }
  s = 7.0/2.0;
  for(int i=1;i<nTerms;++i){
    term3 += a7_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  //my linear parameterization of Chi term using fit to values in Table 1
  double intercept = 1.60986e-6, slope = -1.83031e-9;
  term4 = (intercept + slope * T) * Hi*1000;
  double mag = spont_m + term1 + term2 + term3 + term4;
  return mag;
}

void comboplot(){
  //density of Ni at room temperature in g/cm^3
  const double RHO_Ni = 8.902;
  //magnetic saturation induction of Ni in gaus at H_int=1.2 T
  const double H_sat = 55.24*(4*3.1415927*RHO_Ni)/1000.0;
  
  if(0){
    TGraph *grt = new TGraph();
    for(int i=1;i<100;++i)grt->SetPoint(i-1,i,magnetization(i,286.4));
    grt->Draw("acp");
    return;
  }

  TGaxis::SetMaxDigits(4);
  gStyle->SetPadRightMargin(0.075);
  gStyle->SetPadLeftMargin(0.11);

  TCanvas *c= new TCanvas("c","c", 0, 0, 660, 1200);
  TPad *pad1 = new TPad("pad1","",0,0.5,1,1);//c->Divide(1,2);
  TPad *pad2 = new TPad("pad2","",0,0,1,0.5);
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  double lower_y = 53.7;
  double upper_y = 55.7;
  double lower_limit =  0;
  double upper_limit =  28.3;
  int style[7] = {34,21,8,4,33,34,26};
  int color[9] = {kBlue,kBlack,kGreen+3,kRed,kBlue,kRed,kViolet,kOrange+7,1};
  //int color[7] = {kBlue,kBlue+3,kBlue-4,kBlue-7,kAzure+7,kBlue-5,kViolet+9};
  const double pi = 3.1415927;

  //Danan 1959 
  const int nDan = 53, nDanC = 96;
  double xDanan[nDan] = {0.03762, 0.04104, 0.04831, 0.04155, 0.05200, 0.06006, 0.05769, 0.06668, 0.09876, 0.10033, 0.09959, 0.11865, 0.11786, 0.11508, 0.14805, 0.17440, 0.18943, 0.21865, 0.23027, 0.23542, 0.26829, 0.27045, 0.29991, 0.33670, 0.36304, 0.36515, 0.40488, 0.45738, 0.47139, 0.50980, 0.52030, 0.55084, 0.55351, 0.57810, 0.60518, 0.65078, 0.71251, 0.74352, 0.78085, 0.80354, 0.84345, 0.89143, 0.94504, 0.99506, 1.00468, 1.29168, 1.47658, 1.71303, 1.80852, 1.11157, 0.91709, 0.45777, 0.05533};//internal field H_i

  double yDanan[nDan] = {55.36077, 55.48553, 55.42849, 55.34289, 55.47288, 55.34933, 55.37966, 55.44273, 55.35950, 55.32035, 55.28209, 55.35679, 55.31166, 55.27118, 55.27974, 55.30702, 55.25163, 55.20173, 55.27733, 55.18732, 55.17665, 55.15884, 55.21912, 55.23078, 55.12037, 55.13585, 55.11940, 55.18901, 55.10741, 55.16969, 55.12747, 55.09093, 55.11401, 55.12013, 55.13048, 55.10783, 55.10124, 55.11097, 55.03800, 55.07656, 55.07099, 55.09098, 55.05772, 55.01405, 55.06060, 54.98401, 54.96047, 54.86009, 54.82799, 55.01149, 55.01812, 55.17043, 55.45524};
  double xDananC[nDanC]={0.03755, 0.04318, 0.04640, 0.05528, 0.05851, 0.06821, 0.07387, 0.08196, 0.08924, 0.12485, 0.13214, 0.16614, 0.17314, 0.21128, 0.23179, 0.25125, 0.27124, 0.29486, 0.31358, 0.33483, 0.35481, 0.37480, 0.39479, 0.41478, 0.43476, 0.45475, 0.47474, 0.49473, 0.51396, 0.53653, 0.55652, 0.61480, 0.63405, 0.65401, 0.67646, 0.69721, 0.73644, 0.75643, 0.77642, 0.79316, 0.82352, 0.86183, 0.88212, 0.90055, 0.92179, 0.96177, 0.98176, 1.00301, 1.02174, 1.04173, 1.06172, 1.08170, 1.12587, 1.14712, 1.16711, 1.18710, 1.20708, 1.22707, 1.24706, 1.26705, 1.28679, 1.30702, 1.32701, 1.34700, 1.36699, 1.38697, 1.40696, 1.42695, 1.44694, 1.46692, 1.48691, 1.50690, 1.52870, 1.54869, 1.56868, 1.58866, 1.60865, 1.62864, 1.64862, 1.66861, 1.68860, 1.70859, 1.72857, 1.74856, 1.76854, 1.78853, 1.80852, 1.82850, 1.84849, 1.86847, 1.88846, 1.90845, 1.92843, 1.94842, 1.96840, 1.98566};//internal field H_i
  double yDananC[nDanC] = {55.50374, 55.46171, 55.44393, 55.41240, 55.40189, 55.37764, 55.36632, 55.35177, 55.33964, 55.29677, 55.29111, 55.26603, 55.26157, 55.23914, 55.22892, 55.21967, 55.21083, 55.20218, 55.19471, 55.18540, 55.17871, 55.17227, 55.16682, 55.16062, 55.15525, 55.15078, 55.14596, 55.14061, 55.13547, 55.13043, 55.12710, 55.11970, 55.11671, 55.11324, 55.10998, 55.10609, 55.10170, 55.09676, 55.09161, 55.08776, 55.08310, 55.07447, 55.06987, 55.06543, 55.06034, 55.05174, 55.04546, 55.04031, 55.03497, 55.02943, 55.02455, 55.01918, 55.00607, 55.00034, 54.99554, 54.98976, 54.98332, 54.97663, 54.97117, 54.96514, 54.96057, 54.95324, 54.94688, 54.94085, 54.93457, 54.92821, 54.92202, 54.91566, 54.90971, 54.90401, 54.89797, 54.89138, 54.88377, 54.87749, 54.87105, 54.86427, 54.85758, 54.85048, 54.84387, 54.83726, 54.83090, 54.82397, 54.81686, 54.80910, 54.80183, 54.79383, 54.78582, 54.77847, 54.77029, 54.76311, 54.75502, 54.74792, 54.74007, 54.73231, 54.72397, 54.71741};
  int n = 0;
  for(int i=0;i<nDan;++i)
    xDanan[i] = 1.000/xDanan[i];
  for(int i=0;i<nDanC;++i)
    xDananC[i] = 1.000/xDananC[i];

  TGraph *grDanan = new TGraph(nDan, xDanan, yDanan);
  grDanan->SetMarkerColor(color[n]);
  grDanan->SetMarkerStyle(style[n]);
  TGraph *grDananC = new TGraph(nDanC, xDananC, yDananC);
  grDananC->SetLineColor(color[n]);
  grDanan->SetMarkerColor(color[n]);
  grDanan->SetMarkerStyle(style[n]);
  grDanan->GetXaxis()->SetLimits(lower_limit,upper_limit);
  grDanan->GetYaxis()->SetLimits(213.4,upper_y);
  grDanan->GetYaxis()->SetRangeUser(lower_y,upper_y);
  grDanan->Draw("ap");
  grDanan->SetTitle(Form("Magnetization of Nickel vs Internal Field H_{int}"));
  grDanan->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grDanan->GetXaxis()->SetTitle("H_{int} (kOe)");
  grDanan->GetXaxis()->SetTitleSize(0.042);
  grDanan->GetYaxis()->SetTitleSize(0.042);
  //grDananC->Draw("samec");
  ++n;
  const int nAra = 12;
  double xAraj[nAra] = {0.15827,0.13022,0.24642,0.33257,0.50286,0.68517,0.97567,1.30023,1.74699,2.07155,2.57041,3.31368};//internal field H_i

  double yAraj[nAra] = {55.26579, 55.31196, 55.19958, 55.16550, 55.10737, 55.06632, 54.98316, 54.89801, 54.74568, 54.61536, 54.43293, 53.93944};
  for(int i=0;i<nAra;++i)
    xAraj[i] = 1.000/xAraj[i];
  TGraph *grAraj = new TGraph(nAra, xAraj, yAraj);
  grAraj->SetMarkerColor(color[n]);
  grAraj->SetMarkerStyle(style[n]);
  grAraj->SetLineColor(color[n]);
  grAraj->Draw("samep");
  ++n;

  const int nCra = 8;
  const double demagCrangle = 5.76;
  double xCrangle[nCra] = {4000,5000,6000,7000,8000,9000,10000,10500};//applied field H
  double yCrangle[nCra] = {55.15, 55.075,55.095,55.105,55.095,55.11,55.17,55.105};
  for(int i=0;i<nCra;++i){
    xCrangle[i] -= yCrangle[i]*4*pi*RHO_Ni/demagCrangle;
    xCrangle[i] *= 0.001;
  }
  TGraph *grCrangle = new TGraph(nCra, xCrangle, yCrangle);
  grCrangle->SetMarkerColor(color[n]);
  grCrangle->SetMarkerStyle(style[n]);
  grCrangle->SetMarkerSize(1);
  grCrangle->SetLineColor(color[n]);
  grCrangle->Draw("samep");
  ++n;


  const int nShu = 13;
  const double demagShull = 3.0;
  double xShull[nShu] = {2965.67633, 3479.92475, 3933.06598, 4946.08977, 6003.69102, 6994.61478, 7973.00915, 8983.39121, 10013.08060,19988.04469,29962.78233,39974.32318,49959.31078};//applied field H
  double yShull[nShu] = {54.84565,54.90577, 54.92201, 54.95896, 54.95082, 55.01134, 55.00442, 55.02715, 55.01130, 55.11516, 55.21780, 55.23352, 55.24882};

  for(int i=0;i<nShu;++i){
    xShull[i] -= yShull[i]*4*pi*RHO_Ni/demagShull;
    xShull[i] *= 0.001;
  }
  TGraph *grShull = new TGraph(nShu, xShull, yShull);
  grShull->SetMarkerColor(color[n]);
  grShull->SetMarkerStyle(style[n]);
  grShull->SetMarkerSize(1);
  grShull->SetLineColor(color[n]);
  grShull->Draw("samep");



  TLegend *leg = new TLegend(0.56, 0.13, 0.92, 0.6);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(grDanan, "Henri Danan (1959)","p");
  leg->AddEntry(grAraj, "Araj #it{et al.} (1967)","p");
  leg->AddEntry(grCrangle, "Crangle #it{et al.} (1970)","p");
  leg->AddEntry(grShull, "Shull #it{et al.} (NIST) (2000)","p");
  leg->Draw();
  double x_ = 5, y_ = 54.7, x_err = 0, y_err = y_*0.002;
  TGraphErrors *gr02 = new TGraphErrors(1,&x_,&y_,&x_err,&y_err);
  gr02->SetLineWidth(3);
  gr02->Draw("samep");
  TPaveText *pte = new TPaveText(0.17,0.423,0.33,0.477,"ndc");
  pte->SetBorderSize(0);
  pte->SetShadowColor(0);
  pte->SetTextColor(kBlack);
  pte->SetFillColor(0);
  pte->AddText("0.4%");
  pte->Draw();
  TPaveText *tp1 = new TPaveText(0.2,0.1,0.5,0.25,"ndc");
  tp1->SetShadowColor(0);
  tp1->SetFillColor(0);
  tp1->AddText("Temperature Range");
  tp1->AddText("288_{ }K to 298_{ }K");
  tp1->Draw();
  pad2->cd();
  n=0;
  //Now correct to T=294 K
  double T = 294;
  for(int i=0;i<nDan;++i){//from 288
    double corrDanan = magnetization(xDanan[i], T);
    corrDanan -= magnetization(xDanan[i], 288);
    yDanan[i] += corrDanan;
  }
  for(int i=0;i<nDanC;++i){//from 288
    double corrDanan = magnetization(xDananC[i], T);
    corrDanan -= magnetization(xDananC[i], 288);
    yDananC[i] += corrDanan;
  }
  for(int i=0;i<nAra;++i){//from 298
    double corrAraj = magnetization(xAraj[i], T);
    corrAraj -= magnetization(xAraj[i], 298);
    yAraj[i] += corrAraj;
  }
  for(int i=0;i<nCra;++i){//from 293
    double corrCrangle = magnetization(xCrangle[i], T);
    corrCrangle -= magnetization(xCrangle[i], 293);
    yCrangle[i]+=corrCrangle;
  }
  for(int i=0;i<nShu;++i){//from 298
    double corrShull = magnetization(xShull[i], T);
    corrShull -= magnetization(xShull[i], 298);
    yShull[i]+=corrShull;
  }

  TGraph *grDanan2 = new TGraph(nDan, xDanan, yDanan);
  grDanan2->SetMarkerColor(color[n]);
  grDanan2->SetMarkerStyle(style[n]);
  TGraph *grDananC2 = new TGraph(nDanC, xDananC, yDananC);
  grDananC2->SetLineColor(color[n]);
  grDanan2->SetMarkerColor(color[n]);
  grDanan2->SetMarkerStyle(style[n]);
  grDanan2->GetXaxis()->SetLimits(lower_limit,upper_limit);
  grDanan2->GetYaxis()->SetRangeUser(lower_y,upper_y);
  grDanan2->Draw("ap");
  grDanan2->SetTitle(Form("Magnetization of Nickel at 294 K versus H_{int}"));
  grDanan2->SetTitle("");
  grDanan2->GetXaxis()->SetTitleSize(0.042);
  grDanan2->GetYaxis()->SetTitleSize(0.042);
  grDanan2->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grDanan2->GetXaxis()->SetTitle("H_{int} (kOe)");
  //  grDananC2->Draw("samec");
  ++n;
  TGraph *grAraj2 = new TGraph(nAra, xAraj, yAraj);
  grAraj2->SetMarkerColor(color[n]);
  grAraj2->SetMarkerStyle(style[n]);
  grAraj2->SetLineColor(color[n]);
  grAraj2->Draw("samep");
  ++n;
  TGraph *grCrangle2 = new TGraph(nCra, xCrangle, yCrangle);
  grCrangle2->SetMarkerColor(color[n]);
  grCrangle2->SetMarkerSize(1);
  grCrangle2->SetMarkerStyle(style[n]);
  grCrangle2->SetLineColor(color[n]);
  grCrangle2->Draw("samep");
  ++n;
  TGraph *grShull2 = new TGraph(nShu, xShull, yShull);
  grShull2->SetMarkerColor(color[n]);
  grShull2->SetMarkerStyle(style[n]);
  grShull2->SetLineColor(color[n]);
  grShull2->Draw("samep");
  leg->Draw();
  gr02->Draw("samep");

  TPaveText *pte2 = new TPaveText(0.172,0.523,0.33,0.575,"ndc");
  pte2->SetBorderSize(0);
  pte2->SetShadowColor(0);
  pte2->SetFillColor(0);
  pte2->AddText("0.4%");
  pte2->Draw();
  TPaveText *tp2 = new TPaveText(0.2,0.25,0.5,0.35,"ndc");
  tp2->SetShadowColor(0);
  tp2->SetFillColor(0);
  tp2->AddText("Corrected to 294_{ }K");
  tp2->Draw();
  // TF1 *f = new TF1("f","[0]+[1]*sqrt(x)+[2]*x",100,20000);
  // f->SetParameters(218, 1, 0);
  // grDanan->Fit(f,"r");
  c->SaveAs(Form("NiMagnetization_vs_Hint.pdf"));
  c->SaveAs(Form("../nim/figures/NiMagnetization_vs_Hint.pdf"));

  TCanvas *c2 = new TCanvas("c2","c2",0,600,1000,660);
  c2->SetGrid();
  TGraph *grAll = new TGraph();
  n=0;
  for(int i=0;i<nDan;++i)grAll->SetPoint(i+n, xDanan[i], yDanan[i]);
  n+=nDan;
  for(int i=0;i<nAra;++i)grAll->SetPoint(i+n, xAraj[i], yAraj[i]);
  n+=nAra;
  for(int i=0;i<nCra;++i)grAll->SetPoint(i+n, xCrangle[i], yCrangle[i]);
  n+=nCra;
  for(int i=0;i<nShu;++i)grAll->SetPoint(i+n, xShull[i], yShull[i]);
  n+=nShu;
  grAll->SetMarkerColor(kGray+2);
  grAll->SetMarkerStyle(8);
  grAll->Draw("ap");
  grAll->SetTitle(Form("Magnetization of Nickel at 294 K vs H_{int}"));
  grAll->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grAll->GetXaxis()->SetTitle("H_{int} (kOe)");
  //  grAll->GetYaxis()->SetTitleOffset(1.4);
  grAll->GetXaxis()->SetRangeUser(lower_limit, upper_limit);
  grAll->GetYaxis()->SetRangeUser(lower_y,upper_y);
  gPad->Update();
  TF1 *f = new TF1("f","magnetization(x+[1],294,[0])",0,29);
  f->SetLineColor(kGray);
  f->SetParameters(58,0);
  grAll->Fit(f,"r");
  gStyle->SetOptFit(0);
  printf("Results from fitting all points at once.\n");
  printf("Fit Msat: %f, offset: %f\n",f->GetParameter(0),f->GetParameter(1));
  printf("Evaluated at 14kOe: %f\n",f->Eval(14));
  if(1){
    TF1 *f = new TF1("f","magnetization(x,294,[0])",1.3,20);
    TF1 *f1 = new TF1("f1","magnetization(x+[1],294,[0])",1.0,20);
    TF1 *f2 = new TF1("f2","magnetization(x+[1],294,[0])",1.5,20);
    TF1 *f3 = new TF1("f3","magnetization(x+[1],294,[0])",0.5,15);
    f->SetParameter(0,58);
    f1->SetParameters(58, 0);
    f2->SetParameters(58, 0);
    f3->SetParameters(58, 0);
    f->SetLineColor(kBlue);
    f1->SetLineColor(kBlack);
    f2->SetLineColor(kRed);
    f3->SetLineColor(kGreen+2);
  
    grAll->Fit(f,"r");
    grAll->Fit(f1,"r");
    grAll->Fit(f2,"r");
    grAll->Fit(f3,"r");
    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    f->Draw("same");
    c2->ForceUpdate();
    c2->SaveAs("NiParameterization_vs_Hint.pdf");
    TCanvas *cx = new TCanvas("cx","cx",0,0,660,1200);
    TPad *pad1x = new TPad("pad1x","",0,0.5,1,1);
    TPad *pad2x = new TPad("pad2x","",0,0,1,0.5);
    pad1x->SetBottomMargin(0);
    pad2x->SetTopMargin(0);
    pad1x->Draw();
    pad2x->Draw();

    pad1x->cd();

    TF1 *fx = new TF1("fx","magnetization(x,294,[0])+[1]/x/x",0,30);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *grf[4];
    double param[2] = {0,0};
    grf[0]=(TGraph*)grDanan2->Clone();
    grf[1]=(TGraph*)grAraj2->Clone();
    grf[2]=(TGraph*)grCrangle2->Clone();
    grf[3]=(TGraph*)grShull2->Clone();
    for(int i=0;i<4;++i){
      fx->SetRange(0,30);
      fx->SetParameters(58.4,-0.3);
      fx->SetLineColor(color[i]);
      if(i==2){
	fx->FixParameter(1,0);
	grf[i]->Fit(fx,"rB");
	param[0]+=fx->GetParameter(0)/4.;
	fx->ReleaseParameter(1);
	//	fx->ReleaseParameter(2);
      }else{
	fx->SetParLimits(1,-10,0);
	grf[i]->Fit(fx,"rB");
	param[0]+=fx->GetParameter(0)/4.;
	param[1]+=fx->GetParameter(1)/3.;
      }
      mg->Add(grf[i]);
      fx->Draw("same");
    }
    fx->SetParameters(param[0],param[1]);
    fx->SetLineWidth(3);
    fx->SetLineStyle(10);
    fx->SetLineColor(kBlack);
    leg->AddEntry(fx,"Average","l");
    mg->SetTitle(Form("Magnetization of Nickel at 294 K vs H_{int}"));
    mg->Draw("ap");
    mg->GetYaxis()->SetTitleSize(0.042);
    mg->GetYaxis()->SetTitle("Magnetization (emu/g)");
    mg->GetXaxis()->SetTitle("H_{int} (kOe)");
    mg->GetXaxis()->SetLimits(lower_limit, upper_limit);
    mg->GetXaxis()->SetRangeUser(lower_limit, upper_limit);
    mg->GetYaxis()->SetRangeUser(lower_y, upper_y);
    mg->Draw("ap");
    fx->Draw("same");
    leg->Draw();
    gr02->Draw("samep");
    pte->SetX1NDC(0.19);
    pte->Draw();
    gPad->Update();
    
    pad2x->cd();
    double x[100],xe[100],y[100],ye[100];
    for(int i=0;i<100;++i){
      x[i] = (i+1)*0.2815;
      xe[i] = 0;
      y[i] = fx->Eval(x[i]);
      ye[i] = 0.002*y[i];
    }
    TGraph *grAll2 = new TGraph();
    grAll2->SetMarkerStyle(8);
    double xt, yt;
    for(int i=0;i<grAll->GetN();++i){
      grAll->GetPoint(i, xt, yt);
      grAll2->SetPoint(i, xt, yt);
    }
    TGraphErrors *gr = new TGraphErrors(100,x,y,xe,ye);
    TGraph *gr1 = new TGraphErrors(100,x,y);
    gr1->SetLineColor(kBlue);
    gr1->SetLineWidth(2);
    gr->SetFillColor(kCyan);
    gr->GetXaxis()->SetTitleSize(0.042);
    gr->GetYaxis()->SetTitleSize(0.042);
    gr->SetTitle(Form(""));
    gr->GetYaxis()->SetTitle("Magnetization (emu/g)");
    gr->GetXaxis()->SetTitle("H_{int} (kOe)");
    gr->GetXaxis()->SetLimits(lower_limit, upper_limit);
    gr->GetXaxis()->SetRangeUser(lower_limit, upper_limit);
    gr->GetYaxis()->SetRangeUser(lower_y, upper_y);
    gr->Draw("3A");
    gr1->Draw("samec");
    grAll2->Draw("samep");
    gr02->Draw("samep");
    pte2->Draw();
    TArrow *ar4 = new TArrow(6,53.93,20,53.93,0.025,"<|>");
    ar4->SetAngle(40);
    ar4->SetLineWidth(2);
    ar4->SetFillColor(1);
    ar4->Draw();
    gPad->Update();
    TPaveText *tp3 = new TPaveText(0.34,0.17,0.63,0.235,"ndc");
    tp3->SetBorderSize(0);
    tp3->SetShadowColor(0);
    tp3->SetFillColor(0);
    tp3->AddText("Region of interest");
    tp3->Draw();
    cx->SaveAs("../nim/figures/NiParameterizationErrorBand_vs_Hint.pdf");

    TCanvas *c4 = new TCanvas("c4","c4",0,0,700,500);
    TGraph *grxx = new TGraph();
    fx->SetParameters(param[0],param[1]);
    for(int i=0;i<100;++i){
      grxx->SetPoint(i,6+i*0.140,fx->Eval(6+i*0.140));
    }
    grxx->SetMarkerStyle(8);
    grxx->Draw("ap");
    TF1 *fp2 = new TF1("fp2","pol2",6,20);
    grxx->Fit(fp2);
    double p_0 = fp2->GetParameter(0);
    double p_1 = fp2->GetParameter(1);
    double p_2 = fp2->GetParameter(2);
    TPaveText *pt = new TPaveText(0.4,0.2,0.89,0.35,"ndc");
    pt->AddText("With H_{int} in gaus:");
    pt->AddText(Form("%0.3f%+0.4eH_{int}%+0.4eH_{int}^{2}", p_0, p_1, p_2));
    pt->AddText("With H in tesla:");
    pt->AddText(Form("%0.3f%+0.6fH%+0.6fH^{2}",
		     p_0 - p_1*H_sat + p_2*H_sat*H_sat,
		     (p_1 - 2*p_2*H_sat)*10., p_2*100.));
    printf("%0.3f%+0.4eH_{int}%+0.4eH_{int}^{2}\n", p_0, p_1, p_2);
    printf("%0.3f%+0.6fH%+0.6fH^{2}\n", p_0 - p_1*H_sat + p_2*H_sat*H_sat,
	   (p_1 - 2*p_2*H_sat)*10., p_2*100.);
    printf("Pol2 approx H=2T (Hi=%fgaus) on thin foil: %f (emu/g)\n",
	   20-H_sat, fp2->Eval(20-H_sat));
    pt->SetFillColor(0);
    pt->Draw("same");
    printf("%f  %f\n",fx->Eval(20-H_sat),fp2->Eval(20-H_sat));
  }  
  return;
}
