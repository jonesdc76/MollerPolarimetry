#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "TStyle.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include <sstream>
#include <vector>
const int nVAR = 6;

int Regress(){
  gStyle->SetTitleW(0.95);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.12);
.q  std::vector<std::vector<double>> X(nVAR);
  for(int i=0;i<nVAR;++i)X.push_back(std::vector<double>());
  std::vector<double>Y,x,Yreg, diff,t;
  //ifstream f("Weekend1.55TTempHallProbeLIA032921.txt");
  ifstream f("Overnight1.55TTempHallProbeLIA032521.txt");
  if(!f){cout<<"Exiting\n"; return -1;}
  string line = "";
  int col[nVAR] = {0,4,6,10,12,13};
  cout<<"hi\n";
  getline(f,line);
  int n=0;
  while(getline(f,line)){
    //cout<<line<<"\n";
    std::stringstream ss(line);
    std::string str;
    int var = 0;
    for(int i=0;i<=col[nVAR-1];++i){
      std::getline(ss,str,',');
      if(i==0){
	t.push_back(atof(str.data()));
	X[var].push_back(1);
	++var;
      }
      if(i==2){
	Y.push_back(atof(str.data()));
	cout<<n<<" Y: "<<Y.back()<<" X: ";
      }else if(i==col[var]){
	X[var].push_back(atof(str.data()));
	cout<<" "<<X[var].back();
	++var;
      }
    }
    cout<<endl;
    x.push_back(n);
    ++n;
  }
 std:vector<double>m, v;
  for(int i=0;i<nVAR;++i){
    for(int j=0;j<nVAR;++j){
      double sum = 0, sum2 = 0;
      for(int k=0;k<(int)X[i].size();++k){
	sum += X[i][k]*X[j][k];
	if(i==0)
	  sum2 += X[j][k]*Y[k];
      }
      cout<<sum<<endl;
      m.push_back(sum);
      v.push_back(sum2);
    }
  }
  TMatrixD *XXp = new TMatrixD(nVAR, nVAR, &m[0]);
  TMatrixD *XXpInv = new TMatrixD(nVAR, nVAR, &m[0]);
  TVectorD *XpY = new TVectorD(nVAR, &v[0]);
  XXpInv->Invert();
  XXp->Print();
  XpY->Print();
  TVectorD slopes = (*XXpInv)*(*XpY);
  cout<<"The vector of slopes is\n";
  TGraph *gr = new TGraph((int)Y.size(),&(t[0]),&(Y[0]));
  for(int k=0;k<(int)X[0].size();++k){
    double sum=0;
    for(int i=0;i<nVAR;++i){
      sum += X[i][k]*slopes[i];
    }
    Yreg.push_back(sum);
    diff.push_back((sum-Y.at(k))/Y.at(k)*100);
  }
  TCanvas *c = new TCanvas("c","c",0,0,1200,600);
  c->Divide(2,1);
  c->cd(1);
  TGraph *gr1 = new TGraph((int)Yreg.size(),&(t[0]),&(Yreg[0]));
  TGraph *gr2 = new TGraph((int)diff.size(),&(t[0]),&(diff[0]));
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.1);
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(0.1);
  gr->SetMarkerColor(kRed);
  gr1->SetMarkerColor(kBlue);
  gr->Draw("ap");
  gr->SetTitle("Measured 1F vs Time (Red) and 5-Variable Regression Model of 1F (Blue)");
  gr->GetXaxis()->SetTitle("Time (s)");
  gr->GetYaxis()->SetTitle("1F (V)");
  gr1->Draw("samep");
  c->cd(2);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(0.1);
  gr2->Draw("ap");
  gr2->SetTitle("Percent Difference Between Regression Model and Measured 1F");
  gr2->GetXaxis()->SetTitle("Time (s)");
  gr2->GetYaxis()->SetTitle("Difference (%)");
  slopes.Print();
  f.close();
  return 0;
}
