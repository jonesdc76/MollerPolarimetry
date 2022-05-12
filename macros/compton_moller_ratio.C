#include <TGraphErrors.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
int plot(){ 
  double rat[10] = {0.998999,
		    0.995567,
		    1.002895,
		    0.998401,
		    0.996803,
		    1.001004,
		    0.998372,
		    1.003047,
		    0.996145,
		    1.000115};
  double err[10] = {0.003412489,
		    0.003426534,
		    0.002107380,
		    0.002013224,
		    0.002754979,
		    0.002530698,
		    0.002839242,
		    0.002143168,
		    0.002287486,
    		    0.005245116};
  const char* da[10] = {"02/08/2020","02/24/2020","03/18/2020","08/19/2020","09/04/2020","09/16/2020"};
  double x[10] = {0.95,1.05,1.95,2.05,3,3.95,4.05,4.95,5.05,6}, xe[10]={0,0,0,0,0,0,0,0,0,0};
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TGraphErrors *gr = new TGraphErrors(10,x,rat,xe,err);
  gr->SetMarkerStyle(8);
  gr->Draw("ap");
  gr->GetXaxis()->SetNdivisions(8);
  gr->GetYaxis()->SetTitleSize(0.044);
  gr->GetYaxis()->SetTitleOffset(1.1);
  gr->GetYaxis()->SetLimits(0.99,1.011);
  gr->GetYaxis()->SetRangeUser(0.9905,1.009);
  gr->SetTitle("");
  gr->GetYaxis()->SetTitle("Ratio (Compton/M#lower[-0.1]{#slash{#lower[0.1]{#kern[0.19]{o}}}}ller)");
  //gr->GetYaxis()->SetTitle("Ratio (Compton/M#slash{o}ller)");
  TF1 *f = new TF1("f","pol0",0,1);
  f->SetLineStyle(10);
  f->SetLineWidth(3);
  f->SetLineColor(kBlack);
  gr->Fit(f);
  TPaveText *ptwo = new TPaveText(0.1,0.05,0.95,0.1391,"ndc");
  ptwo->SetFillColor(0);
  ptwo->SetBorderSize(0);
  ptwo->SetShadowColor(0);
  ptwo->Draw();
  TPaveText *pt[6];
  for (int i=0; i<6; i++){
    pt[i] = new TPaveText(0.06+0.1395*(i+1)-0.11,0.013,0.06+0.1395*(i+1),0.1388,"ndc");
    pt[i]->SetFillColor(0);
    pt[i]->SetBorderSize(0);
    pt[i]->SetShadowColor(0);
    TText *tt = pt[i]->AddText(da[i]);
    tt->SetTextAngle(30);
    tt->Draw();
    pt[i]->Draw();
  }
  TPaveText *pto = new TPaveText(0.6,0.8,0.89,0.94,"ndc");
  pto->AddText(Form("#chi^{2}/NDF: %0.1f/%i",f->GetChisquare(),f->GetNDF()));
  pto->AddText(Form("Ratio: %0.5f #pm %0.5f",f->GetParameter(0),f->GetParError(0)));
  pto->SetFillColorAlpha(0,0);
  pto->SetShadowColor(0);
  //  pto->SetTextColor(kRed);
  pto->SetLineColor(0);
  pto->Draw();
  gPad->Update();
  c->SaveAs("CREX_Compton_to_Moller_ratio.pdf");
  cout<<f->GetProb()<<endl;
  return 0;
}
