#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
using namespace std;

int plot(bool is_prex = true, bool by_group = false){
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
  TGraphErrors *grIn= new TGraphErrors(), *grOut = new TGraphErrors();
  int nin = 0, nout = 0, n = 0;
  TPaveText *pt[20];
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
    cout<<date<<endl;
    date += " 00:00:00";
    TDatime da(00,00,00,00,00,00);
    da.Set(date.data());
    hwp = words[2];
    double polfac = (is_prex ? 1615.9827 : 1656.3281);
    bool new_point = (hwp != prevhwp||date != prev_date||by_group)  && !first;
    if(new_point && asym!=0){
      if(TString(prevhwp.data()).Contains("IN")){
	asym_in.push_back(asym/tmpAe);
	asym_in_e.push_back(1/sqrt(tmpAe));
	printf("%i In %0.8f+/-%0.8f\n",prev_group,asym_in.back(),asym_in_e.back());
	double x = double(by_group ? prev_group : da.Convert());
	grIn->SetPoint(nin,x,abs(asym_in.back())*polfac);
	grIn->SetPointError(nin,0,asym_in_e.back()*polfac);
	++nin;
      }else{
	asym_out.push_back(asym/tmpAe);
	asym_out_e.push_back(1/sqrt(tmpAe));
	printf("%i Out %0.8f+/-%0.8f\n",prev_group, asym_out.back(),asym_out_e.back());
	double x = double(by_group ? prev_group : da.Convert());
	grOut->SetPoint(nout,x,abs(asym_out.back()*polfac));
	grOut->SetPointError(nout,0,asym_out_e.back()*polfac);
	++nout;
      }
      // pt[n] = new TPaveText(0.12+n*0.1,0.12,0.12+n*0.1+0.1,0.16,"ndc");
      // pt[n]->SetFillColor(0);
      // pt[n]->SetShadowColor(0);
      // pt[n]->SetBorderSize(0);
      // pt[n]->AddText(Form("%s/%s/%s",month.data(),day.data(),year.data()));
      ++n;
      asym = 0;
      tmpAe = 0;
    }
    if(is_prex && group>2000)continue;
    else if(!is_prex && (group<3000 || group > 4000) )continue;
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
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  c->SetGrid();
  TMultiGraph *mg = new TMultiGraph();
  grIn->SetMarkerColor(kBlue);
  grIn->SetLineColor(kBlue);
  grIn->SetMarkerStyle(8);
  TF1 *f = new TF1("f","pol0",0,1);
  f->SetLineWidth(3);
  f->SetLineStyle(10);
  f->SetLineColor(kRed);
  grOut->Fit(f);
  TPaveText *pto = new TPaveText(0.6,0.78,0.89,0.898,"ndc");
  pto->AddText(Form("#chi^{2}/NDF: %0.1f/%i",f->GetChisquare(),f->GetNDF()));
  pto->AddText(Form("Polarization: %0.2f #pm %0.2f%%",f->GetParameter(0),f->GetParError(0)));
  f->SetLineColor(kBlue);
  grIn->Fit(f);
  TPaveText *pti = new TPaveText(0.6,0.655,0.89,0.7685,"ndc");
  pti->AddText(Form("#chi^{2}/NDF:      %0.1f/%i",f->GetChisquare(),f->GetNDF()));
  pti->AddText(Form("Polarization: %0.2f #pm %0.2f%%",f->GetParameter(0),f->GetParError(0)));
  grOut->SetMarkerColor(kRed);
  grOut->SetLineColor(kRed);
  grOut->SetMarkerStyle(8);
  grOut->SetMarkerSize(2);
  grIn->SetMarkerSize(2);
  mg->Add(grIn);
  mg->Add(grOut);
  mg->Draw("ap");
  if(is_prex)
    mg->SetTitle("M#philler Polarimetry Measurements for PREX-2");
  else
    mg->SetTitle("M#philler Polarimetry Measurements for CREX");
  if(by_group)
    mg->GetXaxis()->SetTitle("Group No.");
  else     mg->GetXaxis()->SetTitle("Measurement No.");
  mg->GetYaxis()->SetTitle("Polarization (%)");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetLimits(88,92.495);
  mg->GetYaxis()->SetRangeUser(88,92.495);
  mg->Draw("ap");
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetNdivisions(503);
  mg->GetXaxis()->SetTimeFormat("%m/%d/%Y");
  mg->GetXaxis()->SetTimeOffset(0,"gmt");
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
  TLegend *lg = new TLegend(0.4,0.75,0.6,0.89);
  lg->SetFillColorAlpha(0,0);
  lg->SetShadowColor(0);
  lg->SetBorderSize(0);
  //  lg->SetTextColor(kRed);
  lg->AddEntry(grOut,"HWP Out","lp");
  // lg->SetTextColor(kBlue);
  lg->AddEntry(grIn,"HWP In","lp");
  lg->Draw();
  //for(int i=0;i<n;++i)pt[i]->Draw();
  return 0;
}
