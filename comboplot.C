/*Plots the saturation magnetization data for Fe from numerous papers including

Sanford et al 1941 "A determination of the magnetic saturation induction of iron at room temperature"

Danan 1959 "On the interpretation of the magnetization measurements of pure polycrystalline iron and nickel in the vicinity of saturation"

Araj et al 1967 "Electrical Resistivity Studies of Chromium-Rich Chromium-Cobalt Alloys

Crangle et al 1971 "The Magnetization of Pure Iron and Nickel"

Behrendt et al 1972 "Saturation Magnetization of Polycrystalline Iron"

Plots are converted where possible to be versus applied field H

 */
#include <iostream>
#include <cstdio>
#include "TGaxis.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"

//Returns theoretical magnetization value for given Hi in Oersted, T in Kelvin
//Curves from Pauthenet Mar 1982 "Spin Waves in nickel, iron and yttrium-iron
//garnet" equation 9 and Table 1
//Adjustable parameter is spontaneous magnetization (nominally 222.678 emu) 
double magnetization( double Hi, double T, double spont_m = 222.678){
  //calculate 307e-6*T^(3/2)*F(3/2, 1.378*Hi/T) from equation (9) at T and Hi
  const double a3_2 = 307e-6, a5_2 = 22.8e-8, b = 1.378e-4;
  double s = 3.0/2.0, term1 = 0, term2 = 0, term3 = 0;
  for(int i=1;i<10000;++i){
    term1 += a3_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  s = 5.0/2.0;
  for(int i=1;i<10000;++i){
    term2 += a5_2*pow(T,s)*pow(i,-s)*exp(-i*b*Hi/T);
  }
  //my linear parameterization of Chi term using fit to values in Table 1
  double intercept = 3.644e-6, slope = 5.0434e-10;
  term3 = (intercept + slope * T) * Hi;
  double mag = spont_m - term1 - term2 + term3;
  return mag;
}

void comboplot(bool use_Hi = true){ 

  if(0){
    TGraph *grt = new TGraph();
    for(int i=1;i<100;++i)grt->SetPoint(i-1,i*1000,magnetization(i*1000,286.4));
    grt->Draw("acp");
    return;
  }

  TGaxis::SetMaxDigits(4);
  gStyle->SetPadRightMargin(0.075);
  gStyle->SetPadLeftMargin(0.11);

  TCanvas *c= new TCanvas("c","c",0,0,660,1200);
  c->Divide(1,2);
  c->cd(1);
  double upper_y = 219;
  double lower_limit = (use_Hi ? 0:2000);
  double upper_limit = (use_Hi ? 20000:40000);
  int style[7] = {34,21,8,4,33,34,26};
  int color[9] = {kGray+2,kBlack,kGreen+3,kRed,kBlue,kRed,kViolet,kOrange+7,1};
  //int color[7] = {kBlue,kBlue+3,kBlue-4,kBlue-7,kAzure+7,kBlue-5,kViolet+9};
  const double pi = 3.1415926;

  //Weiss and Forrer 1929
  const int nWeis = 3;
  int n = 0;
  const double demagWeiss = 7.24;//taken to be that of a prolate ellipsoid of 
                                 // ratio 9:4 (see Weiss and Forrer pg 300)
  double xWeiss[nWeis] = {1e4,1e4,1e4};//applied field H

  //Y-value uses eqs 1 and 30 from Weiss and Forrer with a=2.6 to evaluate at 
  //H=10000 Oe. Same value given 3 times for extra weighting.
  double yWeiss[nWeis] = {217.76*(1-2.6/10000.), 217.76*(1-2.6/10000.),
			  217.76*(1-2.6/10000.)};

  if(use_Hi)
    for(int i=0;i<nWeis;++i)xWeiss[i] -= yWeiss[i]*4*pi*7.874/demagWeiss;
  TGraph *grWeiss = new TGraph(nWeis, xWeiss, yWeiss);
  grWeiss->SetMarkerColor(color[n]);
  grWeiss->SetMarkerSize(1.5);
  grWeiss->SetMarkerStyle(style[n]);
  ++n;

  //Sanford 1941 single point
  const int nSan = 3;
  const double demagSanford = 32.1;//taken to be that of a cylinder of diameter
                                   //0.6 mm and length 8 mm (see Sanford pg 7 
                                   // and Table II and eq 2b  "Simple and 
                                   //approximate expressions of demagnetizing 
                                   //factors of uniformly magnetized 
                                   //rectangular rod and cylinder" M. Sato)

  double ySanford[nSan] = {217.99,217.99,217.99}, 
    xSanford[nSan] = {10000, 10000, 10000};
  if(use_Hi)
    for(int i=0;i<nSan;++i)  
      xSanford[i] -= ySanford[i]*4*pi*7.874/demagSanford;
  TGraph *grSanford = new TGraph(nSan, xSanford, ySanford);
  grSanford->SetMarkerColor(color[n]);
  grSanford->SetMarkerSize(1.1);
  grSanford->SetMarkerStyle(style[n]);
  ++n;

  //Danan 1959 (Demagnetizing factor same as Weiss and Forrer)
  const int nDan = 33, nDanC = 89;
  const double demagDanan = demagWeiss;
  double xDanan[nDan] = {0.03551, 0.04017, 0.09126, 0.04759, 0.05587, 0.06866, 0.08866, 0.12895, 0.16360, 0.19827, 0.21007, 0.20912, 0.25013, 0.22837, 0.29310, 0.36874, 0.3696480378516221, 0.4224972238925199, 0.5801450159543885, 0.6540048367912149, 0.7086403876455537, 0.7424452569161288, 0.7825550722665766, 0.8208960792145261, 0.8673240079530542, 0.9184241365537642, 1.0085609448860369, 1.0752138590307982, 1.1964847676189554, 1.2932205919143949, 1.4018157718017707, 1.5981908118523331,  1.7526586989821666};
  double yDanan[nDan] = {218.00568, 217.91716, 217.85066, 217.80878, 217.74197, 217.71134, 217.73859, 217.56358, 217.52407, 217.47372, 217.50634, 217.54248, 217.51386, 217.44141, 217.38405, 217.34483,217.34484082918930, 217.32353612451007, 217.25057936156114, 217.17701690243464, 217.18464077388660, 217.09633092956807, 217.06589798847423, 216.99207876683770, 217.0086794506599, 216.9277235479751,  216.95006188634858, 216.8113870881376, 216.69479386579002, 216.5328425584958,  216.35471192124194, 215.88986643960376, 215.43556138818542};//internal field H_i

  double xDananC[nDanC]={  0.02825, 0.04113, 0.05672, 0.0901963, 0.1102606, 0.1259643, 0.1448999, 0.1641936, 0.1875700, 0.2346587, 0.2545963, 0.2742439, 0.2919695, 0.3118779, 0.3325178, 0.3525502, 0.3725799, 0.3926154, 0.4126468, 0.4392051, 0.4599958, 0.4800294, 0.4999142, 0.5200964, 0.5401308, 0.5601636, 0.5801950, 0.6002320, 0.6202680, 0.6382160, 0.7004492, 0.6803710, 0.7243529, 0.7690656, 0.8078771, 0.8265169, 0.8443119, 0.8643436, 0.8821105, 0.9026114, 0.9178829, 0.9323183, 0.9536185, 0.9736579, 0.9944564, 1.0137445, 1.0337787, 1.0545762, 1.0920832, 1.1121254, 1.1321690, 1.1522111, 1.1722569, 1.1926044, 1.2123472, 1.2323942, 1.2524413, 1.2724889, 1.2925421, 1.3126005, 1.3326529, 1.3527078, 1.3727594, 1.3928138, 1.4128678, 1.4329222, 1.4529747, 1.4730294, 1.4930842, 1.5131387, 1.5331954, 1.5532513, 1.5733073, 1.5940976, 1.6134158, 1.6306757, 1.6517158, 1.6717763, 1.6918344, 1.7118955, 1.7319554, 1.7511069, 1.7702589, 1.7903196, 1.8103815, 1.8304434, 1.8495940, 1.8678364, 1.8869925};//internal field H_i
  double yDananC[nDanC] = {  217.9821, 217.8810, 217.7871, 217.6644, 217.6244, 217.5956, 217.5693, 217.5394, 217.5190, 217.4762, 217.4620, 217.4440, 217.4319, 217.4203, 217.4057, 217.3919, 217.3803, 217.3640, 217.3509, 217.3325, 217.3187, 217.3039, 217.2867, 217.2745, 217.2590, 217.2449, 217.2319, 217.2143, 217.1975, 217.1846, 217.1349, 217.1514, 217.1116, 217.0741, 217.0412, 217.0216, 217.0060, 216.9928, 216.9715, 216.9466, 216.9333, 216.9191, 216.8986, 216.8791, 216.8587, 216.8335, 216.8183, 216.7988, 216.7547, 216.7329, 216.7099, 216.6880, 216.6632, 216.6396, 216.6146, 216.5888, 216.5629, 216.5366, 216.5057, 216.4704, 216.4402, 216.4078, 216.3783, 216.3463, 216.3147, 216.2829, 216.2525, 216.2203, 216.1881, 216.1561, 216.1222, 216.0891, 216.0559, 216.0173, 215.9923, 215.9564, 215.9196, 215.8826, 215.8477, 215.8102, 215.7738, 215.7362, 215.6982, 215.6611, 215.6230, 215.5849, 215.5480, 215.5099, 215.4685};
  
  for(int i=0;i<nDan;++i)
    xDanan[i] = 1000/xDanan[i]+(use_Hi ? 0:4*pi*yDanan[i]*7.874/demagDanan);
  for(int i=0;i<nDanC;++i)
    xDananC[i] = 1000/xDananC[i]+(use_Hi ? 0:4*pi*yDananC[i]*7.874/demagDanan);

  TGraph *grDanan = new TGraph(nDan, xDanan, yDanan);
  grDanan->SetMarkerColor(color[n]);
  grDanan->SetMarkerStyle(style[n]);
  TGraph *grDananC = new TGraph(nDanC, xDananC, yDananC);
  grDananC->SetLineColor(color[n]);
  grDanan->SetMarkerColor(color[n]);
  grDanan->SetMarkerStyle(style[n]);
  grDanan->GetXaxis()->SetLimits(lower_limit,upper_limit);
  grDanan->GetYaxis()->SetLimits(213.4,upper_y);
  grDanan->GetYaxis()->SetTitleOffset(1.5);
  grDanan->GetYaxis()->SetRangeUser(213.4,upper_y);
  grDanan->Draw("ap");
  grDanan->SetTitle(Form("Magnetization of Iron vs %s",(use_Hi ? "Internal Field H_{int}":"Applied Field H")));
  grDanan->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grDanan->GetXaxis()->SetTitle((use_Hi ? "H_{int} (Oe)" : "H (Oe)"));
  //grDananC->Draw("samec");
  grSanford->Draw("samep");
  ++n;
  const int nAra = 15;
  const double demagAraj = 3.0;
  double xAraj[nAra] = {1.83963, 1.59581, 1.38579, 1.24863, 1.10280, 0.92821, 0.71221, 0.61922, 0.51771, 0.47053, 0.41331, 0.33462, 0.28456, 0.25172, 0.23176};//internal field H_i

  double yAraj[nAra] = {215.06297, 215.87245, 216.27994, 216.54179, 216.75649, 216.88718, 217.12381, 217.21786, 217.36050, 217.42667, 217.48839, 217.56186, 217.62213, 217.69717, 217.77076};

  for(int i=0;i<nAra;++i)
    xAraj[i] = 1000/xAraj[i] + (use_Hi ? 0:4*pi*yAraj[i]*7.874/demagAraj);
  TGraph *grAraj = new TGraph(nAra, xAraj, yAraj);
  grAraj->SetMarkerColor(color[n]);
  grAraj->SetMarkerStyle(style[n]);
  grAraj->SetLineColor(color[n]);
  grAraj->Draw("samepc");
  ++n;

  const int nCra = 8;
  const double demagCrangle = 5.76;
  double xCrangle[nCra] = {4000,5000,6000,7000,8000,9000,10000,10500};//applied field H
  double yCrangle[nCra] = {213.5, 217.5, 217.65, 217.7, 217.55, 217.75, 217.6, 217.7};
  if(use_Hi)
    for(int i=0;i<nCra;++i)
      xCrangle[i] -= yCrangle[i]*4*pi*7.874/demagCrangle;
  TGraph *grCrangle = new TGraph(nCra, xCrangle, yCrangle);
  grCrangle->SetMarkerColor(color[n]);
  grCrangle->SetMarkerStyle(style[n]);
  grCrangle->SetMarkerSize(1.5);
  grCrangle->SetLineColor(color[n]);
  grCrangle->Draw("samepl");
  grWeiss->Draw("samep");
  ++n;

  const int nNasa = 26;
  const double demagNasa = 3.0;
  double xNasa[nNasa] = {0.18380,0.28653,0.37976,0.50229,0.60599,0.69708,0.80355,0.90037,0.99829,1.10086,1.19899,1.29509,1.39033,1.49674,1.60303,1.69727,1.79041,1.89866,2.00191,2.09468,2.19793,2.29533,2.40148,2.50363,2.59619,2.70248};
  double yNasa[nNasa] = {217.38795, 217.61164, 217.70140, 217.73676, 217.78903, 217.80291, 217.86079, 217.90370, 217.92225, 217.99606, 218.03522, 218.00788, 218.02925, 218.34755, 218.38763, 218.31064, 218.25894, 218.49106, 218.49929, 218.41200, 218.42024, 218.52029, 218.67934, 218.44683, 218.60404,218.64412};

  for(int i=0;i<nNasa;++i)
    xNasa[i] = 10000*xNasa[i] + (use_Hi ? 0:4*pi*yNasa[i]*7.874/demagNasa);
  TGraph *grNasa = new TGraph(nNasa, xNasa, yNasa);
  grNasa->SetMarkerColor(color[n]);
  grCrangle->SetMarkerSize(1.5);
  grNasa->SetMarkerStyle(style[n]);
  grNasa->SetLineColor(color[n]);
  grNasa->Draw("samep");
  ++n;
  TLegend *leg = new TLegend(0.6, 0.18, 0.9, 0.63);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(grWeiss, "Weiss #it{et al.} (1929)","p");
  leg->AddEntry(grSanford, "NIST (1941)","p");
  leg->AddEntry(grDanan, "Henri Danan (1959)","p");
  leg->AddEntry(grAraj, "Araj #it{et al.} (1967)","p");
  leg->AddEntry(grCrangle, "Crangle #it{et al.} (1970)","p");
  leg->AddEntry(grNasa, "NASA(1972)","p");
  leg->Draw();
  c->cd(2);
  n=0;
  //Now correct to T=294 K
  double corrWeiss = -0.142;//from 288
  for(int i=0;i<nWeis;++i)yWeiss[i]+=corrWeiss;
  double corrSanford = 0.096;//from 298
  for(int i=0;i<nSan;++i)ySanford[i]+=corrSanford;
  double corrDanan = -0.02;//from 293
  for(int i=0;i<nDan;++i)yDanan[i]+=corrDanan;
  for(int i=0;i<nDanC;++i)yDananC[i]+=corrDanan;
  double corrAraj = 0.096;//from 298
  for(int i=0;i<nAra;++i)yAraj[i]+=corrAraj;
  double corrCrangle = -0.02;//from 293
  for(int i=0;i<nCra;++i)yCrangle[i]+=corrCrangle;
  double corrNasa = 0.12;//from 299
  for(int i=0;i<nNasa;++i)yNasa[i]+=corrNasa;
  TGraph *grWeiss2 = new TGraph(nWeis, xWeiss, yWeiss);
  grWeiss2->SetMarkerColor(color[n]);
  grWeiss2->SetMarkerSize(1.5);
  grWeiss2->SetMarkerStyle(style[n]);
  ++n;
  TGraph *grSanford2 = new TGraph(nSan, xSanford, ySanford);
  grSanford2->SetMarkerColor(color[n]);
  grSanford2->SetMarkerSize(1);
  grSanford2->SetMarkerStyle(style[n]);
  ++n;

  TGraph *grDanan2 = new TGraph(nDan, xDanan, yDanan);
  grDanan2->SetMarkerColor(color[n]);
  grDanan2->SetMarkerStyle(style[n]);
  TGraph *grDananC2 = new TGraph(nDanC, xDananC, yDananC);
  grDananC2->SetLineColor(color[n]);
  grDanan2->SetMarkerColor(color[n]);
  grDanan2->SetMarkerStyle(style[n]);
  grDanan2->GetXaxis()->SetLimits(lower_limit,upper_limit);
  grDanan2->GetYaxis()->SetLimits(213.4,upper_y);
  grDanan2->GetYaxis()->SetTitleOffset(1.5);
  grDanan2->GetYaxis()->SetRangeUser(213.4,upper_y);
  grDanan2->Draw("ap");
  grDanan2->SetTitle(Form("Magnetization of Iron at 294 K versus %s",(use_Hi ? "H_{int}":"H")));
  grDanan2->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grDanan2->GetXaxis()->SetTitle((use_Hi ? "H_{int} (Oe)":"H (Oe)"));
  //  grDananC2->Draw("samec");
  grSanford2->Draw("samep");
  ++n;
  TGraph *grAraj2 = new TGraph(nAra, xAraj, yAraj);
  grAraj2->SetMarkerColor(color[n]);
  grAraj2->SetMarkerStyle(style[n]);
  grAraj2->SetLineColor(color[n]);
  grAraj2->Draw("samep");
  ++n;
  TGraph *grCrangle2 = new TGraph(nCra, xCrangle, yCrangle);
  grCrangle2->SetMarkerColor(color[n]);
  grCrangle2->SetMarkerSize(1.5);
  grCrangle2->SetMarkerStyle(style[n]);
  grCrangle2->SetLineColor(color[n]);
  grCrangle2->Draw("samepl");
  ++n;
  TGraph *grNasa2 = new TGraph(nNasa, xNasa, yNasa);
  grNasa2->SetMarkerColor(color[n]);
  grNasa2->SetMarkerStyle(style[n]);
  grNasa2->SetLineColor(color[n]);
  grNasa2->Draw("samep");
  grWeiss2->Draw("samep");
  leg->Draw();
  // TF1 *f = new TF1("f","[0]+[1]*sqrt(x)+[2]*x",100,20000);
  // f->SetParameters(218, 1, 0);
  // grDanan->Fit(f,"r");
  c->SaveAs(Form("Magnetization_vs_%s.pdf",(use_Hi ? "Hint":"H")));
  if(0) return;
  TCanvas *c2 = new TCanvas("c2","c2",0,600,1000,660);
  c2->SetGrid();
  TGraph *grAll = new TGraph();
  n=0;
  for(int i=0;i<nWeis;++i)grAll->SetPoint(i, xWeiss[i], yWeiss[i]);
  n+=nWeis;
  for(int i=0;i<nSan;++i){
    cout<<xSanford[i]<<" "<<ySanford[i]<<endl;
    grAll->SetPoint(i+n, xSanford[i], ySanford[i]);
  }
  n+=nSan;
  for(int i=0;i<nDan;++i)grAll->SetPoint(i+n, xDanan[i], yDanan[i]);
  n+=nDan;
  for(int i=0;i<nAra;++i)grAll->SetPoint(i+n, xAraj[i], yAraj[i]);
  n+=nAra;
  for(int i=0;i<nCra;++i)grAll->SetPoint(i+n, xCrangle[i], yCrangle[i]);
  n+=nCra;
  for(int i=0;i<nNasa;++i)grAll->SetPoint(i+n, xNasa[i], yNasa[i]);
  n+=nNasa;
  grAll->SetMarkerColor(kGray+2);
  grAll->SetMarkerStyle(8);
  grAll->Draw("ap");
  grAll->SetTitle(Form("Magnetization of Iron at 294 K vs %s",(use_Hi ? "H_{int}":"H")));
  grAll->GetYaxis()->SetTitle("Magnetization (emu/g)");
  grAll->GetXaxis()->SetTitle((use_Hi ? "H_{int} (Oe)":"H (Oe)"));
  grAll->GetXaxis()->SetTitleOffset(1.4);
  grAll->GetXaxis()->SetRangeUser(500,20000);
  grAll->GetYaxis()->SetRangeUser(215,218.5);
  gPad->Update();
  gStyle->SetOptFit(0);


  // TF1 *f = new TF1("f","magnetization(x,294,[0])",1000,20000);
  // TF1 *f1 = new TF1("f1","magnetization(x+[1],294,[0])",1000,20000);
  // TF1 *f2 = new TF1("f2","magnetization(x+[1],294,[0])",500,20000);
  // TF1 *f3 = new TF1("f2","magnetization(x+[1],294,[0])",500,15000);
  // f->SetParameter(0,222);
  // f1->SetParameters(222, 2130);
  // f2->SetParameters(222, 1931);
  // f3->SetParameters(222, 1925);
  // f->SetLineColor(kBlue);
  // f1->SetLineColor(kBlack);
  // f2->SetLineColor(kRed);
  // f3->SetLineColor(kGreen+2);
  
  //grAll->Fit(f,"r");
  // grAll->Fit(f1,"r");
  // grAll->Fit(f2,"r");
  // grAll->Fit(f3,"r");
  // f1->Draw("same");
  // f2->Draw("same");
  // f3->Draw("same");
  // f->Draw("same");
  c2->ForceUpdate();
  // if(use_Hi)
  //   c2->SaveAs("CombinedFit_vs_Hint.pdf");
  // else
  //   c2->SaveAs("CombinedFit_vs_H.pdf");
 
  return;
}
