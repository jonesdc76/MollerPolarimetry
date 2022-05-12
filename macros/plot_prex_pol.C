#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
using namespace std;

int plot(bool by_group = false){
  gStyle->SetOptFit(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin(0.075);
  ifstream file("PCrexAsym.dat");
  string line, date, hwp, tmp, prevhwp=" ", prev_date=" ";
  double asym=0, asym_e=0;
  int group, prev_group = 0;
  string delim = ",";
  vector<string>words;
  vector<double>asym_in, asym_out, asym_in_e, asym_out_e;
  double x[100], xe[100], tmpAe=0;
  getline(file, line);
  bool first = true;
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TGraphErrors *grIn= new TGraphErrors(), *grOut = new TGraphErrors();
  int nin = 0, nout = 0, n = 0;
  double xmin = 0, xmax = 0;
  while(getline(file, line)){
    size_t pos=0;
    words.clear();
    string::iterator end_pos = remove(line.begin(),line.end(),' ');
    line.erase(end_pos, line.end());
    while((pos=line.find(delim, 0)) != string::npos){
      words.push_back(line.substr(0,pos));
      //cout<<line.substr(0,pos)<<"    "<<line.length()<<endl;
      line.erase(0, pos+delim.length());
    }
    pos = line.length();
    words.push_back(line.substr(0,pos));
    group = atoi(words[0].data());
    date = words[1];
    hwp = words[2];
    double polfac = 1615.9827; //1656.3281);
    bool new_point = (hwp != prevhwp||date != prev_date|| by_group);
    if(new_point && asym != 0){
      if(TString(prevhwp.data()).Contains("IN")){
	asym_in.push_back(asym/tmpAe);
	asym_in_e.push_back(1/sqrt(tmpAe));
	printf("%i In %0.8f+/-%0.8f\n",prev_group,asym_in.back(),asym_in_e.back());
	double x = double(by_group ? prev_group : nin+1);
	cout<<x<<endl;
	grIn->SetPoint(nin,x,abs(asym_in.back())*polfac);
	grIn->SetPointError(nin,0,asym_in_e.back()*polfac);
	++nin;
      }else{
	asym_out.push_back(asym/tmpAe);
	asym_out_e.push_back(1/sqrt(tmpAe));
	printf("%i Out %0.8f+/-%0.8f\n",prev_group, asym_out.back(),asym_out_e.back());
	double x = double(by_group ? prev_group : nout+1);
	grOut->SetPoint(nout,x,abs(asym_out.back()*polfac));
	grOut->SetPointError(nout,0,asym_out_e.back()*polfac);
	++nout;
      }
      ++n;
      asym = 0;
      tmpAe = 0;
    }
    if(group>2000)continue;
    double scale = (group>1045 ? 1.0 : 1.0110);
    asym_e = atof(words[4].data());
    asym_e *= scale;
    asym += atof(words[3].data())*scale/pow(asym_e, 2);
    tmpAe += 1/pow(asym_e, 2);
    prev_group = group;
    prevhwp = hwp;
    prev_date = date;
    first = false;
  }
  const char* da[8] = {"08/04/2019","08/10/2019","08/18/2019","08/21/2019","08/26/2019","08/31/2019","09/04/2019","09/08/2019"};
  TPaveText *pt[8];
  if(!by_group){
    for (int i=0; i<8; i++){
      pt[i] = new TPaveText(0.07+(0.1015*i),0.01,0.07+0.1015*(i+1),0.114,"ndc");
      pt[i]->SetFillColor(0);
      pt[i]->SetBorderSize(0);
      pt[i]->SetShadowColor(0);
      TText *tt = pt[i]->AddText(da[i]);
      tt->SetTextAngle(30);
      tt->Draw();
    }
  }

  TMultiGraph *mg = new TMultiGraph();
  grIn->SetMarkerColor(kBlue);
  grIn->SetLineColor(kBlue);
  grIn->SetMarkerStyle(8);
  TF1 *f = new TF1("f","pol0",0,1);
  f->SetLineWidth(3);
  f->SetLineStyle(10);
  f->SetLineColor(kRed);
  grOut->Fit(f);
  TPaveText *pto = new TPaveText(0.6,0.83,0.89,0.948,"ndc");
  pto->AddText(Form("#chi^{2}/NDF: %0.1f/%i",f->GetChisquare(),f->GetNDF()));
  pto->AddText(Form("Polarization: %0.2f #pm %0.2f%%",f->GetParameter(0),f->GetParError(0)));
  f->SetLineColor(kBlue);
  grIn->Fit(f);
  TPaveText *pti = new TPaveText(0.6,0.705,0.89,0.8185,"ndc");
  pti->AddText(Form("#chi^{2}/NDF:      %0.1f/%i",f->GetChisquare(),f->GetNDF()));
  pti->AddText(Form("Polarization: %0.2f #pm %0.2f%%",f->GetParameter(0),f->GetParError(0)));
  grOut->SetMarkerColor(kRed);
  grOut->SetLineColor(kRed);
  grOut->SetMarkerStyle(8);
  grOut->SetMarkerSize(1);
  grIn->SetMarkerSize(1);
  mg->Add(grIn);
  mg->Add(grOut);
  mg->SetTitle("");
  mg->Draw("ap");
  if(by_group)
    mg->GetXaxis()->SetTitle("Group No.");
  else     mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("Polarization (%)");
  mg->GetXaxis()->SetNdivisions(15);
  if(!by_group){
    mg->GetXaxis()->SetLimits(0,11);
    mg->GetXaxis()->SetRangeUser(0.5,8.5);
  }
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleOffset(1.05);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetLimits(88,92.495);
  mg->GetYaxis()->SetRangeUser(88.2,92.1);
  mg->Draw("ap");
  pto->SetFillColorAlpha(0,0);
  pto->SetShadowColor(0);
  pto->SetTextColor(kRed);
  pto->SetLineColor(0);
  pti->SetFillColorAlpha(0,0);
  pti->SetShadowColor(0);
  pti->SetTextColor(kBlue);
  pti->SetLineColor(0);
  pti->Draw();
  pto->Draw();
  
  TLegend *lg = new TLegend(0.15,0.77,0.35,0.91);
  lg->SetFillColorAlpha(0,0);
  lg->SetShadowColor(0);
  lg->SetBorderSize(0);
  lg->AddEntry(grOut,"HWP Out","lp");
  lg->AddEntry(grIn,"HWP In","lp");
  lg->Draw();
  if(!by_group) for (int i=0; i<8; i++)pt[i]->Draw();
  gPad->Update();
  c->SaveAs((by_group ? "PREXpol_vs_group.pdf":"PREXpol_vs_date.pdf"));
  return 0;
}
