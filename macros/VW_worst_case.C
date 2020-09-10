#include <iostream>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
struct fit{
  double off;
  double amp;
  double ang;
};
const double pi = 3.141592653590;
double theta = pi/3.0;//angle of birefringent optic w.r.t. the x-axis
double Delta = 0.02;//retardance in radians

double Re_x1(double theta, double Delta){
  return (pow(cos(theta),2)+pow(sin(theta),2)*cos(Delta)+cos(theta)*sin(theta)*sin(Delta))/sqrt(2);
}
double Im_x1(double theta, double Delta){
  return (pow(sin(theta),2)*sin(Delta)+cos(theta)*sin(theta)-cos(theta)*sin(theta)*cos(Delta))/sqrt(2);
}
double Re_x2(double theta, double Delta){
  return (cos(theta)*sin(theta)-cos(theta)*sin(theta)*cos(Delta)-pow(cos(theta),2)*sin(Delta))/sqrt(2);
}
double Im_x2(double theta, double Delta){
  return (-cos(theta)*sin(theta)*sin(Delta)+pow(sin(theta),2)+pow(cos(theta),2)*cos(Delta))/sqrt(2);
}
TString func(double theta, double Delta){
  TString str = Form("pow(cos(x)*cos(x)*%f+cos(x)*sin(x)*%f,2)+"
		     "pow(cos(x)*cos(x)*%f+cos(x)*sin(x)*%f,2)+"
		     "pow(cos(x)*sin(x)*%f+sin(x)*sin(x)*%f,2)+"
		     "pow(cos(x)*sin(x)*%f+sin(x)*sin(x)*%f,2)",
		     Re_x1(theta,Delta), Re_x2(theta,Delta), Im_x1(theta,Delta),
		     Im_x2(theta,Delta), Re_x1(theta,Delta), Re_x2(theta,Delta),
		     Im_x1(theta,Delta), Im_x2(theta,Delta));
  return str;
}

//Vacuum window birefringence
void matrix_bif(double theta, double delta, TMatrixDSym* m, bool re){
  double arr_re[4], arr_im[4], arr[4];
    arr_re[0] = pow(cos(theta),2)+pow(sin(theta),2)*cos(delta);
    arr_re[1] = cos(theta)*sin(theta)*(1-cos(delta));
    arr_re[2] = cos(theta)*sin(theta)*(1-cos(delta));
    arr_re[3] = pow(sin(theta),2)+pow(cos(theta),2)*cos(delta);
    arr_im[0] = pow(sin(theta),2)*sin(delta);
    arr_im[1] = -cos(theta)*sin(theta)*sin(delta);
    arr_im[2] = -cos(theta)*sin(theta)*sin(delta);
    arr_im[3] = +pow(cos(theta),2)*sin(delta);
  if(re){
    arr[0] = arr_re[0]*cos(delta/2.0)+arr_im[0]*sin(delta/2.0);
    arr[1] = arr_re[1]*cos(delta/2.0)+arr_im[1]*sin(delta/2.0);
    arr[2] = arr_re[2]*cos(delta/2.0)+arr_im[2]*sin(delta/2.0);
    arr[3] = arr_re[3]*cos(delta/2.0)+arr_im[3]*sin(delta/2.0);
  }
  else{
    arr[0] = arr_im[0]*cos(delta/2.0)-arr_re[0]*sin(delta/2.0);
    arr[1] = arr_im[1]*cos(delta/2.0)-arr_re[1]*sin(delta/2.0);
    arr[2] = arr_im[2]*cos(delta/2.0)-arr_re[2]*sin(delta/2.0);
    arr[3] = arr_im[3]*cos(delta/2.0)-arr_re[3]*sin(delta/2.0);
  }
  m->SetMatrixArray(arr);
}

double normsq(TVectorD* re, TVectorD* im){
  return pow((*re)(0),2)+pow((*re)(1),2)+pow((*im)(0),2)+pow((*im)(1),2);
}

fit birefringence(double retardance = 0, double axis_angle = pi/3){
  gStyle->SetOptFit(0);
  bool verbose = 0;
  const double delta_init = 0.1;//initial Delta on beam before vacuum window
  // const int N=18;
  // TF1 *fTsq[N];
  // for(int i=N-1;i>=0;--i){
  //   fTsq[i] = new TF1(Form("fTsq%i",N),func(theta, 0.1*i).Data(),0,2*pi);
  //   fTsq[i]->SetLineColor(kGreen+i);
  //   fTsq[i]->Draw((i==N-1?"":"same"));
  // }
  double param[3] = {0.5,1,0.2};;
  TMatrixDSym M_re = TMatrixDSym(2);
  TMatrixDSym M_im = TMatrixDSym(2);

  TVectorD V_re_init = TVectorD(2);
  V_re_init(0) = 1/sqrt(2);
  V_re_init(1) = -sin(delta_init)/sqrt(2);
  TVectorD V_im_init = TVectorD(2);
  V_im_init(0) = 0;
  V_im_init(1) = cos(delta_init)/sqrt(2);
  TGraph *gr;
  TF1 *fcos = new TF1("fcos","[0]+[1]*cos(2*(x-[2]))",0,2*pi);
  fcos->FixParameter(0,0.5);
  fcos->SetParameter(1,1);
  fcos->SetParameter(2,0.2);
  double ret = retardance;
  matrix_bif(axis_angle,ret,&M_re,1);
  matrix_bif(axis_angle,ret,&M_im,0);
  TVectorD V_re = M_re*V_re_init - M_im*V_im_init;
  TVectorD V_im = M_re*V_im_init + M_im*V_re_init;
  
  if(verbose){
    printf("Re\n[%+f  %+f]\n[%+f   %+f]\n\n",
	   M_re(0,0),M_re(0,1),M_re(1,0),M_re(1,1));
    printf("[Im\n%+f  %+f]\n[%+f   %+f]\n",
	   M_im(0,0),M_im(0,1),M_im(1,0),M_im(1,1));
    printf("Re:[%+f]  Im: [%+f]\n", V_re(0), V_im(0));
    printf("   [%+f]      [%+f]\n", V_re(1), V_im(1));
  }

  //Analyze with linear polarizer
  TMatrixDSym Mpol = TMatrixDSym(2);
  const int n = 500;
  double x[n], y[n];
  for(int j=0;j<n;++j){
    x[j] = j*2*pi/double(n-1);
    Mpol(0,0) = pow(cos(x[j]),2);
    Mpol(0,1) = cos(x[j])*sin(x[j]);
    Mpol(1,0) = cos(x[j])*sin(x[j]);
    Mpol(1,1) = pow(sin(x[j]),2);
    TVectorD vre = Mpol*V_re;
    TVectorD vim = Mpol*V_im;
    y[j] = normsq(&vre, &vim);
  }

  gr = new TGraph(n,x,y);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.2);
  gr->SetMarkerColor(kBlue);   
  gr->SetLineColor(kBlue);   
  gr->SetLineWidth(1);
  fcos->SetLineColor(kRed);   
  gr->Fit(fcos,"q");
  // if(fcos->GetParameter(1)<0){
  //   fcos->SetParameter(1, fabs(fcos->GetParameter(1))+0.01);
  //   fcos->SetParameter(2, fcos->GetParameter(2)+pi);
  //   gr->Fit(fcos,"q");
  // }
  param[0] = fcos->GetParameter(0);
  param[1] = fcos->GetParameter(1);
  param[2] = fcos->GetParameter(2);
  //  printf("Ret: %f, Offset: %f, Amp: %f, Angle offset: %f\n",ret,fcos->GetParameter(0), fcos->GetParameter(1), fcos->GetParameter(2));

  //gr->Draw("alp");
  fit xx;
  xx.off = fcos->GetParameter(0);
  xx.amp = fabs(fcos->GetParameter(1));
  xx.ang = fcos->GetParameter(2);
  return xx;
}

int VW_worst_case(double ret = 0.1){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleW(1);
  const int N=7;
  double col[N] = {kViolet,kViolet+10,kAzure+10,kCyan+2,kGreen+3, kGreen};
  double style[N] = {8,21,22,23,29,33, 36};
  double retard[N] = {-0.2,-0.1,0.1,0.2,0.3,0.4};
  TMultiGraph *mg = new TMultiGraph();
  TGraph *grp[N];
  TLegend *tl = new TLegend(0.84,0.5,0.99,0.89);
  tl->SetHeader("#Delta=-0.1 for All","C");
  for(int j=0;j<N-1;++j){
    grp[j] = new TGraph();
    grp[j]->SetMarkerStyle(style[j]);
    grp[j]->SetMarkerSize(1.3);
    grp[j]->SetLineStyle(j+1);
    grp[j]->SetLineWidth(2);
    grp[j]->SetLineColor(col[j]);
    grp[j]->SetMarkerColor(col[j]);
    //grp[j]->GetXaxis()->SetLimits(-95.,110.);
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,700);
    for(int i=0;i<61;++i){
      fit rslt = birefringence(retard[j],i*3*pi/180.);
      grp[j]->SetPoint(i,i*3-90,rslt.amp/rslt.off);
    }
    tl->AddEntry(grp[j],Form("#Delta_{vw}=%+0.1f",retard[j]),"lp");
    mg->Add(grp[j]);
  }
  TCanvas *c = new TCanvas("c","c",0,0,800,600);
  mg->SetTitle("DOLP vs Angle Between Birefringence Axes of Vacuum Window (#Delta_{vw}) and #Delta from Upstream Optics");
  mg->Draw("alp");
  mg->GetXaxis()->SetLimits(-95.,120.);
  mg->GetXaxis()->SetTitle("Angle (deg)");
  mg->GetYaxis()->SetTitle("DOLP");
  tl->Draw();
  gPad->Update();
  c->SaveAs("../Vacuum_window_delta_vs_angle.png");
  return 0;
}
