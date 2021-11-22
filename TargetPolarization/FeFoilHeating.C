#include "TF1.h"
#include <iostream>
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TString.h"

///////////////////
//Donald C. Jones//
//Nov. 2021      //
///////////////////

////////////////////////////////////////////////////////////////////////////////////////
//FeFoilHeating() calculates and graphs the temperature differnce in a thin circular  //
//Fe foil between its edge held at a fixed temperature T0 and inside a circular       //
//Gaussian-distributed or uniformly rastered electron beam.                           //
//                                                                                    //
//Arguments:                                                                          //
//  beam_cur: beam current in Amperes                                                 //
//  beam_r:   1 sigma beam spot size radius in cm  (beam_r^2=sigma_x^2+sigma_y^2)     //
//  beam_E:   beam energy in GeV                                                      //
//  T0:       ambient (Hall) temperature in Kelvin taken as foil boundary temperature //
//  foil_r:   foil radius in cm (default is 1/2")                                     //
//  uniform:  uniform charge distribution? Otherwise, Gaussian assumed.               //
//                                                                                    //
//Returns:                                                                            //
//the foil temperature difference in degrees K between T0 at the foil edge and the    //
//temperature at the 1-sigma beam radius r_beam.                                      //
//NOTE: it is helpful to recall that for a 2D circular Gaussian distribution the      //
//      since the width is the quadrature sum of the x and y widths this is not the   //
//      same as the 1-sigma width of the projected 1D distribution. The projected     //
//      1D distribution width will be sqrt(2) smaller than the beam_r, the 2D 1-sigma //
//      width. The volume between r=0 and the n-sigma are as follows (where sigma     //
//      is the quadrature sum of sigma_x and sigma_y):                                //
//      1 sigma = 39.35%,  2 sigma = 86.47%, 3 sigma = 98.89%, 4 sigma = 99.97%.      //
////////////////////////////////////////////////////////////////////////////////////////

double FeFoilHeating(double beam_cur = 1e-6, double beam_r=15e-3, double beam_E = 11, double T0 = 294, double foil_r = 0.635, bool uniform = 0){
  gStyle->SetStatY(0.7);
  gStyle->SetStatH(0.2);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleW(0.9);

  bool save_plots = 1;
  const double rho = 7.87;//density of Fe
  const double sigma = 5.670e-12;//Stefan Boltzman constant W/(cm^2 K^4)
  const double Cp = 0.45;//Fe specific heat capacity in J/(g K)
  const double echarge = 1.602e-19;//Coulombs per electron
  const double PI = 3.1415927;//pi obviously

  

  //Use ESTAR data to estimate energy loss as a function of electron energy
  //----------------------------------------------------------------------------------
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  double beam_en[10]={1,2,3,4,5,6,7,8,9,10};//beam energy in GeV
  double stop_en[10]={1.878,1.928,1.957,1.977,1.993, //collision stopping power 
		      2.006,2.017,2.027,2.035,2.043};//in (MeV cm^2/g) using ESTAR
  TGraph *grStop = new TGraph(10,beam_en,stop_en);
  grStop->SetTitle("Electron Stopping Power for Fe vs Beam Energy (ESTAR Data)");
  grStop->SetMarkerStyle(8);
  grStop->Draw("ap");
  grStop->GetXaxis()->SetTitle("Electron Energy");
  grStop->GetYaxis()->SetTitle("Stopping Power (MeV cm^{2}/g)");
  gPad->Update();
  TF1 *fStop = new TF1("fStop","pol5",0,1);//use fit to give continuous function
  grStop->Fit(fStop);  
  double alpha = echarge*fStop->Eval(beam_E)*1e6;//Collision stopping power in (Jcm^2/g)
  cout<<"Stopping power "<<alpha<<" (J cm^2/g)"<<endl;
  if(save_plots)
    c->SaveAs("StoppingPower.pdf");

  

  //Calculate the energy dependent thermal conductivity of Fe using data either from
  //https://www.efunda.com/materials/elements/TC_Table.cfm?Element_ID=Fe
  //or
  //https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
  //----------------------------------------------------------------------------------
  bool data_efunda = 1;
  TCanvas *ct = new TCanvas("ct","ct",0,0,800,600);
  double temp[4] = {250,300,350,400};
  double cond[4] = {0.865,0.802,0.744,0.695};//www.efunda.com
  TGraph *grC = new TGraph(4,temp,cond);
  grC->SetTitle("Fe Thermal Conductivity vs. Temperature");
  grC->SetMarkerStyle(8);
  grC->Draw("ap");
  grC->GetXaxis()->SetTitle("Temperature (k)");
  grC->GetYaxis()->SetTitle("Thermal Conductivity (W/cm K)");
  TF1 *fCond = new TF1("fCond","pol2",0,1);
  grC->Fit(fCond);
  gPad->Update();
  if(!data_efunda)//www.engineeringtoolbox.com
    fCond = new TF1("fCond","0.835-0.001102*(x-273)",0,1);
  double slope = uniform ? 19.5 : 17;
  double guessTemp = T0+slope*beam_cur/1e-6;//starting guess for final foil temperature
  double kappa = fCond->Eval(guessTemp);
  cout<<"Conductivity at "<<guessTemp<<" K is "<<kappa<<endl;
  if(save_plots)
    ct->SaveAs("FeThermalCond.pdf");


  
  //Integral of f(r) gives delta T. Create the integrand f(r) 
  //----------------------------------------------------------------------------------
  double gam = beam_cur/echarge*rho*alpha/kappa/PI/pow(beam_r,2)/(uniform ? 1.0 : 2.0);
  double C = -beam_cur/echarge*alpha*rho/2.0/PI/kappa;
  TF1 *f = new TF1("f",Form("%e/x*exp(-x*x/%e)+%e/x",
			    beam_r*beam_r*gam,2*beam_r*beam_r,C),0,foil_r);


  
  //Improve thermal conductivity estimate using the calculated temperature.
  //Temperature at 1.3*beam_r is a good estimate of the average temperature
  //weighted by a Gaussian beam spot charge distribution. For a uniform distribution
  //0.7*beam_r is a good estimate.
  //-----------------------------------------------------------------------------------
  double r_est = (uniform ? 0.7 : 1.3)*beam_r;
  if(uniform)
    guessTemp = gam*pow(beam_r,2)/2.0*log(foil_r/beam_r)+gam/4*(pow(beam_r,2)-pow(r_est,2))+T0;
  else
    guessTemp = f->Integral(foil_r, r_est)+T0;
  kappa = fCond->Eval(guessTemp);
  gam = beam_cur/echarge*rho*alpha/kappa/PI/pow(beam_r,2)/(uniform ? 1.0 : 2.0);
  C = -beam_cur/echarge*alpha*rho/2.0/PI/kappa;

  cout<<"Conductivity re-calculated at "<<guessTemp<<" K is "<<kappa<<endl;
  f = new TF1("f",Form("%e/x*exp(-x*x/%e)+%e/x",beam_r*beam_r*gam,2*beam_r*beam_r,C),0,foil_r);



  //Graph resulting temperature profile by integrating f(r)dr. Make points red inside
  //beam spot radius (2 sigma if Gaussian).
  //-----------------------------------------------------------------------------------
  const int N=1000;
  double r[N], T[N], dT[N],ri[N],Ti[N], dTi[N];
  int n=0, ni=0;
  double rp = foil_r;
  double red_zone = uniform ? beam_r : 2*beam_r;
  for(int i=0;i<N/2;++i){
    r[i]=rp;
    if(uniform){
      if(rp<red_zone)
	dT[i] = gam*pow(beam_r,2)/2.0*log(foil_r/beam_r)+gam/4.0*(pow(beam_r,2)-pow(rp,2));
      else
	dT[i] = gam*pow(beam_r,2)/2.0*log(foil_r/rp);
    }else{
      dT[i] = f->Integral(foil_r,rp);
    }
    T[i] = dT[i]+T0;
    if(rp<red_zone){
      ri[ni]=rp;
      Ti[ni]=T[i];
      dTi[ni]=dT[i];
      ++ni;
    }
    rp*=0.95;
    ++n;
    if(rp<0.00001)break;
  }
  for(int i=0;i<n;++i){
    r[i+n]=-r[n-i-1];
    dT[i+n] = dT[n-i-1];
    T[i+n] = T[n-i-1];
  }
  for(int i=0;i<ni;++i){
    ri[i+ni]=-ri[ni-i-1];
    dTi[i+ni] = dTi[ni-i-1];
    Ti[i+ni] = Ti[ni-i-1];
  }
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  TGraph *grdT = new TGraph(2*n,r,dT);
  grdT->SetMarkerStyle(8);
  grdT->SetLineWidth(6);
  grdT->SetMarkerSize(0.3);
  grdT->Draw("acp");
  grdT->SetTitle(Form("Foil #DeltaT Profile vs Radial Distance from Foil Center"));
  grdT->GetXaxis()->SetTitle("Radial Distance from Foil Center (cm)");
  grdT->GetYaxis()->SetTitle("#DeltaT (K)");
  TGraph *gridT = new TGraph(2*ni,ri,dTi);
  gridT->SetMarkerStyle(8);
  gridT->SetMarkerColor(kRed);
  gridT->SetLineColor(kRed);
  gridT->SetLineWidth(6);
  gridT->SetMarkerSize(0.4);
  gridT->Draw("samep");
  TPaveText *pt = new TPaveText(0.6,0.3,0.89,0.6,"ndc");
  pt->SetFillColor(0);
  pt->SetShadowColor(0);
  pt->SetBorderSize(0);
  pt->AddText(Form("Beam Energy: %0.1f GeV",beam_E));
  pt->AddText(Form("Beam Current: %0.1f #muA", beam_cur*1e6));
  TString str = Form("Beam Spot Size 1#sigma Radius: %0.1f #mum",beam_r*1e4);
  if(uniform)
    str = Form("Beam Spot Size Radius: %0.1f #mum",beam_r*1e4);
  pt->AddText(str.Data());
  pt->AddText((char*)(uniform ? "Beam Spot Profile: Uniform" : "Beam Spot Profile: Gaussian"));
  pt->AddText(Form("Foil Radius: %0.2f cm",foil_r));
  pt->Draw();
  TLegend *lg = new TLegend(0.62,0.76,0.89,0.89);
  if(uniform){
    lg->AddEntry(grdT,"Outside beam spot","lp");
    lg->AddEntry(gridT,"Inside beam spot","lp");
  }else{
    lg->AddEntry(grdT,"Outside 2#sigma beam spot","lp");
    lg->AddEntry(gridT,"Inside 2#sigma beam spot","lp");
  }
  lg->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",0,0,800,600);
  TGraph *gr = new TGraph(2*n,r,T);
  gr->SetMarkerStyle(8);
  gr->SetLineWidth(6);
  gr->SetMarkerSize(0.3);
  gr->Draw("acp");
  gr->SetTitle(Form("Foil Temperature Profile vs Radial Distance from Foil Center"));
  gr->GetYaxis()->SetTitle("Foil Temperature (K)");
  gr->GetXaxis()->SetTitle("Radial Distance from Foil Center (cm)");
  gr->GetYaxis()->SetRangeUser(T0,T0+grdT->GetYaxis()->GetXmax());
  TGraph *gri = new TGraph(2*ni,ri,Ti);
  gri->SetMarkerStyle(8);
  gri->SetMarkerColor(kRed);
  gri->SetLineColor(kRed);
  gri->SetLineWidth(2);
  gri->SetMarkerSize(0.4);
  gri->Draw("samecp");
  lg->Draw();
  pt->Draw();
  

  //Integrate f(r) weighted by the beam charge distribution to find the average delta T 
  //-----------------------------------------------------------------------------------
  gStyle->SetOptFit(0);
  TF1 *fGaus = new TF1("fGaus","[0]*exp(-x*x/(2*[1]*[1]))+[2]",-2*beam_r,2*beam_r);
  fGaus->SetParameters((guessTemp-T0)/2.,beam_r,T0);cout<<(guessTemp-T0)/2.<<endl;
  fGaus->SetLineWidth(2);
  fGaus->SetLineColor(kRed);
  gr->Fit(fGaus,"r");
  if(uniform){
    fGaus->SetRange(-beam_r,beam_r);
    //fGaus->FixParameter(0,fGaus->GetParameter(0));
    fGaus->FixParameter(2,fGaus->GetParameter(2));
    gr->Fit(fGaus,"r");
  }
  TString fstr = Form("x*exp(-x*x/2./%e)/%e",beam_r*beam_r,beam_r*beam_r);
  if(uniform)fstr = Form("2*x/%e",beam_r*beam_r);
  TString func = Form("(%e*exp(-x*x/(2*%e))+%e)*%s",
		      fGaus->GetParameter(0),pow(fGaus->GetParameter(1),2),
		      fGaus->GetParameter(2), fstr.Data());
  TF1 *fAvgT = new TF1("fAvgT",func.Data(),0,1);
  fAvgT->SetNpx(1000);
   //fAvgT->Draw();
  if(uniform)
    cout<<"dT at 0.7 beam radius is "<<fGaus->Eval(beam_r*0.7)<<endl;
  else
    cout<<"dT at 1.3 sigma is "<<f->Integral(foil_r,beam_r*1.3)<<endl;

  

  //Return average temperature, weighted by the beam spot charge distribution.
  //-----------------------------------------------------------------------------------
  c1->SetGrid();
  c1->cd();
  double avg = fAvgT->Integral(0, (uniform ? 1.0 : 10.0 ) * beam_r);
  TPaveText *pt1 = new TPaveText(0.12,0.74,0.48,0.82,"ndc");
  pt1->SetFillColor(0);
  pt1->SetShadowColor(0);
  //pt1->SetBorderSize(0);
  pt1->SetTextColor(kRed);
  pt1->AddText(Form("<#DeltaT> Charge-weighted over Beam Spot"));
  pt1->AddText(Form("%0.2f K",avg-T0));
  pt1->Draw();
  gPad->Update();
  if(save_plots)
    c1->SaveAs(Form("FoilHeatingdT%s.pdf",(char*)(uniform ? "Uniform":"")));
  c2->SetGrid();
  c2->cd();
  TPaveText *pt2 = new TPaveText(0.12,0.74,0.48,0.82,"ndc");
  pt2->SetFillColor(0);
  pt2->SetShadowColor(0);
  //pt2->SetBorderSize(0);
  pt2->SetTextColor(kRed);
  pt2->AddText(Form("<T> Charge-weighted over Beam Spot"));
  pt2->AddText(Form("%0.2f K",avg));
  pt2->Draw();
  if(save_plots)
    c2->SaveAs(Form("FoilHeatingT%s.pdf",(char*)(uniform ? "Uniform":"")));
  cout<<"Total correction to magnetization for Fe: "<<-0.0238*(avg-T0)<<" emu/g"<<endl;
  return avg;
}
