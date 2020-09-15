#include <iostream>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TCanvas.h>

double delta(double pin, double pout, double diff, double dvw){
  return sqrt(1-pow(pout+dvw,2))-sqrt(1-pow(pin+dvw,2))-diff;
}

int solveDeltavw(){
  gStyle->SetOptStat("rMEn");
  TRandom r(0);
  TH1D *h = new TH1D("h","h",1000,-0.3,0);
  TH1D *hpout = new TH1D("hpout","hpout",1000,0.96,1.0);
  TH1D *hpin = new TH1D("hpin","hpin",1000,0.96,1.0);
  double eps = 0.000001;
  double pout = 0.164, pin = -0.053, diff = 0.0134;
  for(int i=0;i<1e7;++i){
    //    TF1 f("f",Form("sqrt(1-pow(%f+x,2))-sqrt(1-pow(%f+x,2))-%f",pout,pin,r.Gaus(0.0134,0.0014)),0,1);
    double dvw=0, dvw2 = 0.05, val2 = delta(pin,pout,diff,dvw2), dvw1 = 0.0, val1 = delta(pin,pout,diff,dvw1);
    int n=0;
    while(abs(val2)>eps){
      dvw = -val2*(dvw2-dvw1)/(val2 - val1)+dvw2;
      dvw1 = dvw2;
      dvw2 = dvw;
      val1 = val2;
      val2 = delta(pin,pout,diff,dvw);
      ++n;
    }
    if(fabs(pout-0.164)<0.00001 && fabs(pin+0.053)<0.00001 && fabs(diff-0.0134)<0.00001)printf("%f  %f  %f\n",dvw, sqrt(1-pow(pout+dvw,2)), sqrt(1-pow(pin+dvw,2)));
    pout = r.Gaus(0.164,0.02), pin = r.Gaus(-0.053,0.02), diff = r.Gaus(0.0134,0.0014);
    if(i%1000000==0)printf("%i  %i  %f\n",i, n, dvw);
    hpout->Fill(sqrt(1-pow(pout+dvw,2)));
    hpin->Fill(sqrt(1-pow(pin+dvw,2)));
    h->Fill(dvw);
  }
  TCanvas *c = new TCanvas("c","c",0,0,1500,600);
  c->Divide(2,1);
  c->cd(1);
  h->Draw();
  c->cd(2);
  hpin->SetLineColor(kRed);
  hpout->Draw();
  hpin->Draw("sames");
  return 0;
}
