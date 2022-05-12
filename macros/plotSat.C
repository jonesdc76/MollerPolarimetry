#include "TGraph.h"
#include "TMultiGraph.h"
#include <iostream>
#include "TLegend.h"
#include "TPad.h"
#include "TAxis.h"

void norm(TGraph* gr, bool reverse = 0, double tweak = 1){
  double scale = 0;
  double n=5;
  int ni = 0;
  for(int i=reverse? 0:gr->GetN()-1;;--i){
    double x,y;
    gr->GetPoint(i,x,y);
    scale += y/n;
    ++ni;
    cout<<i<<" "<<ni<<" "<<y<<endl;
    if(ni>n-0.00001)break;
    if(reverse)i+=2;
  }
  for(int i=0;i<gr->GetN();++i){
    double x,y;
    gr->GetPoint(i,x,y);
    cout<<i<<" "<<x<<" "<<y<<endl;
    gr->SetPoint(i,x,y/scale*tweak);
  }
  cout<<scale<<endl;
}

int plotSat(){
  TMultiGraph *mg = new TMultiGraph();
  TGraph *grp = new TGraph("NewREpolishedFeFoilRamp0to4.56TRet110030921_An2_136degPEMoffMinFoilAngleNeg1.535","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  norm(grp);
  grp->SetMarkerStyle(8);
  grp->SetMarkerColor(kRed);
  TGraph *grp2 = new TGraph("NewREpolishedFeFoilRamp4.56to1.8TRet110030921_An2_136degPEMoffMinFoilAngleNeg0.535","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  norm(grp2,1,1.00);
  grp2->SetMarkerStyle(8);
  grp->SetMarkerSize(0.5);
  grp2->SetMarkerSize(0.5);
  grp2->SetMarkerColor(kBlue);
  TGraph *grp3 = new TGraph("NewREpolishedFeFoilRamp1.8to4.56TRet110030921_An2_136degPEMoffMinFoilAngleNeg2.535","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  norm(grp3);
  grp3->SetMarkerStyle(8);
  grp3->SetMarkerSize(0.5);
  grp3->SetMarkerColor(kGreen+2);
  TGraph *grp4 = new TGraph("NewREpolishedFeFoilRamp4.56to0TRet110030921_An2_136degPEMoffMinFoilAngleNeg2.535","%*lg, %lg, %*lg, %*lg, %*lg, %*lg, %*lg, %*lg, %lg");
  norm(grp4,1,1.00);
  grp4->SetMarkerStyle(8);
  grp4->SetMarkerSize(0.5);
  grp4->SetMarkerColor(kMagenta);
  mg->Add(grp);
  mg->Add(grp2);
  mg->Add(grp3);
  mg->Add(grp4);
  mg->SetTitle("Saturation Curve on Fe Foil Found in SBU Box with Label 99.99% Purity");
  mg->Draw("ap");
  mg->GetXaxis()->SetRangeUser(1.8,4.8);
  mg->GetYaxis()->SetRangeUser(0.8,1.03);
  mg->GetXaxis()->SetTitle("B-field (T)");
  mg->GetYaxis()->SetTitle("1F/2F");

  gPad->Update();
  TLegend *lg = new TLegend(0.65,0.2,0.89,0.3);
  lg->AddEntry(grp,"-1.535#circ","p");
  lg->AddEntry(grp2,"-0.535#circ","p");
  lg->AddEntry(grp3,"-2.535#circ","p");
  lg->AddEntry(grp4,"-2.535#circ","p");
  lg->Draw();

  return 0;
}
