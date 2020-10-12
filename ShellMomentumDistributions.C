#include <iostream>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
const double Z = 26; //atomic number of element in target
const double alpha = 1/137.0, me = 511.0;
double Nn[4] = {2, 8, Z-12, 2};

void momentum_scale(double *mom){
  for(int n=0;n<4;++n){
    double sum = 0;
    for(int j=0;j<n;++j)sum += Nn[j];
    mom[n] = (Z -0.5*(Nn[n]-1)-sum)*alpha*me;
    cout<<mom[n]<<" ";
  }
      cout<<endl;
}

double get_mom(int n, int l, double p){
  double mom = 0;
  switch (n){
  case 1:
    mom = 10.185916/pow(p*p+1,4);
    break;
  case 2:
    if(l==0)mom = 325.94932*pow(4*p*p-1,2)/pow(4*p*p+1,6);
    else if(l==1)mom = 1738.3964*p*p/pow(4*p*p+1,6);
    else {
      cout<<"Invalid n, l selection\n";
      mom = 0;
    }
    break;
  case 3:
    if(l==0)
      mom = 275.01974*pow(4*pow((9*p*p-1)/(9*p*p+1),2)-1,2)/pow(9*p*p+1,4);
    else if(l==1)
      mom = 79205.686*p*p*pow((9*p*p-1)/pow(9*p*p+1,4),2);
    else if(l==2)
      mom = 570280.94*pow(p,4)/pow(9*p*p+1,8);
    else{
      cout<<"Invalid n, l selection\n";
      mom = 0;
    }
    break;
  case 4:
    //mom = 651.89865*pow(8*pow((16*p*p-1)/(16*p*p+1),2)-4,2)/pow(16*p*p+1,6)*pow(16*p*p-1,2);
    //mom = 651.89865*pow(16*p*p-1,2)*pow(8*pow((16*p*p-1)/(16*p*p+1),2)-4,2)/pow(16*p*p+1,6);
    mom = 651.899 * pow(-1 + 16 * p*p,2) * pow(-4 + 8 * pow(-1 + 16 * p*p,2) / pow(1 + 16 * p*p,2) ,2)/pow(1 + 16 * p*p,6);    break;
  default:
    cout<<"Invalid n, l selection\n";
    mom = 0;
  }
  return mom;
}

int ShellMomentumDistributions(){
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  c->SetLogy();c->SetGrid();
  double mom_scale[4];
  momentum_scale(mom_scale);
  TMultiGraph *mg = new TMultiGraph();
  TGraph *gr1 = new TGraph();
  gr1->SetLineColor(kRed);
  TGraph *gr2 = new TGraph();
  gr2->SetLineColor(kGreen+2);
  TGraph *gr3 = new TGraph();
  gr3->SetLineColor(kBlue);
  TGraph *gr4 = new TGraph();
  gr4->SetLineColor(kBlack);
  TGraph *gr5 = new TGraph();
  gr5->SetLineColor(kBlack);
  TGraph *gr6 = new TGraph();
  gr6->SetLineColor(kGray);
  gr1->SetLineWidth(2);
  gr2->SetLineWidth(2);
  gr3->SetLineWidth(2);
  gr4->SetLineWidth(2);
  double i=0;
  for(int c=0;c<=20000;++c){
    gr1->SetPoint(c,i,get_mom(1,0,i/mom_scale[0])*pow(i/mom_scale[0],2)/mom_scale[0]);
    gr2->SetPoint(c,i,0.25*get_mom(2,0,i/mom_scale[1])*pow(i/mom_scale[1],2)/mom_scale[1]+0.75*get_mom(2,1,i/mom_scale[1])*pow(i/mom_scale[1],2)/mom_scale[1]);
    gr3->SetPoint(c,i,2/Nn[2]*get_mom(3,0,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]+6/Nn[2]*get_mom(3,1,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]+(Nn[2]-8)/Nn[2]*get_mom(3,2,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]);
    gr4->SetPoint(c,double(i),get_mom(4,0,double(i)/mom_scale[3])*pow((double)i/mom_scale[3],2)/mom_scale[3]);
    gr5->SetPoint(c,double(i),(2*get_mom(1,0,i/mom_scale[0])*pow(i/mom_scale[0],2)/mom_scale[0]
		  +2.00*get_mom(2,0,i/mom_scale[1])*pow(i/mom_scale[1],2)/mom_scale[1]
		  +6.00*get_mom(2,1,i/mom_scale[1])*pow(i/mom_scale[1],2)/mom_scale[1]
		  +2.00*get_mom(3,0,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]
		  +6.00*get_mom(3,1,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]
		  +3.79*get_mom(3,2,i/mom_scale[2])*pow(i/mom_scale[2],2)/mom_scale[2]
		  +2.00*get_mom(4,0,i/mom_scale[3])*pow(i/mom_scale[3],2)/mom_scale[3])/23.79);
    gr6->SetPoint(c,i,gr5->Integral(0,c));
    i+=0.01;
    if(i>200)break;
  }
  mg->SetTitle(Form("Electron Momentum Distributions for Z=%i",int(Z)));
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  //  gr10->Draw("al");
  mg->GetYaxis()->SetLimits(1e-5,1);
  mg->Draw("al");
  mg->GetYaxis()->SetTitle("Probability/P_{e} (c/keV)");
  mg->GetXaxis()->SetTitle("P_{e} (keV/c)");
  mg->GetYaxis()->SetLimits(1e-5,1);
  mg->GetYaxis()->SetRangeUser(1e-5,1);
  gPad->Update();
  //gr5->Draw("alp");
  //gr6->Draw("alp");
  return 0;
}
