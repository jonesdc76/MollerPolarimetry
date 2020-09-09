#include <iostream>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <TStyle.h>

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
  double arr[4];
  if(re){
    arr[0] = pow(cos(theta),2)+pow(sin(theta),2)*cos(delta);
    arr[1] = cos(theta)*sin(theta)*(1-cos(delta));
    arr[2] = cos(theta)*sin(theta)*(1-cos(delta));
    arr[3] = pow(sin(theta),2)+pow(cos(theta),2)*cos(delta);
  }
  else{
    arr[0] = pow(sin(theta),2)*sin(delta);
    arr[1] = -cos(theta)*sin(theta)*sin(delta);
    arr[2] = -cos(theta)*sin(theta)*sin(delta);
    arr[3] = pow(cos(theta),2)*sin(delta);
  }
  m->SetMatrixArray(arr);
}

double normsq(TVectorD* re, TVectorD* im){
  return pow((*re)(0),2)+pow((*re)(1),2)+pow((*im)(0),2)+pow((*im)(1),2);
}

double birefringence(double retardance = 0, double axis_angle = pi/3){
  gStyle->SetOptFit(0);
  bool verbose = 0;
  const double delta_init = 0.0;//initial Delta on beam before vacuum window
  TCanvas *c = new TCanvas("c","c",0,0,1700,700);
  c->Divide(2,1);
  c->cd(1);
  const int N=18;
  TF1 *fTsq[N];
  for(int i=N-1;i>=0;--i){
    fTsq[i] = new TF1(Form("fTsq%i",N),func(theta, 0.1*i).Data(),0,2*pi);
    fTsq[i]->SetLineColor(kGreen+i);
    fTsq[i]->Draw((i==N-1?"":"same"));
  }
  c->cd(2);
  TMatrixDSym M_re = TMatrixDSym(2);
  TMatrixDSym M_im = TMatrixDSym(2);

  TVectorD V_re_init = TVectorD(2);
  V_re_init(0) = 1/sqrt(2);
  V_re_init(1) = -sin(delta_init)/sqrt(2);
  TVectorD V_im_init = TVectorD(2);
  V_im_init(0) = 0;
  V_im_init(1) = cos(delta_init)/sqrt(2);
  TMultiGraph *mg = new TMultiGraph();
  const int NN=500;
  const double step = 3.5/double(NN-1);
  TGraph *gr[NN];
  double fit[3][NN], retard[NN], pol[NN];
  TF1 *fcos = new TF1("fcos","[0]+[1]*cos(2*(x-[2]))",0,2*pi);
  fcos->FixParameter(0,0.5);
  fcos->SetParameter(1,1);
  fcos->SetParameter(2,0.2);
  for(int i=NN-1;i>=0;--i){
    double ret = i*step;
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

    gr[i] = new TGraph(n,x,y);
    gr[i]->SetMarkerStyle(8);
    gr[i]->SetMarkerSize(0.2);
    gr[i]->SetMarkerColor(i==0? kRed:kBlack);   
    gr[i]->SetLineColor(i==0? kRed:kRed);   
    gr[i]->SetLineWidth(i==0? 1:1);
    fcos->SetLineColor(i==0||i==NN-1? kRed:kGreen);   
    gr[i]->Fit(fcos,"q");
    if(fcos->GetParameter(1)<0){
      fcos->SetParameter(1,0.25);
      gr[i]->Fit(fcos,"q");
    }
    if(fcos->GetParameter(2)<-4){
      fcos->SetParameter(2,fcos->GetParameter(2)+pi);
      gr[i]->Fit(fcos,"q");
    }
    retard[i]=ret;
    fit[0][i] = fcos->GetParameter(0);
    fit[1][i] = fcos->GetParameter(1);
    fit[2][i] = fcos->GetParameter(2);
    pol[i] = fcos->GetParameter(1)/fcos->GetParameter(0);
    printf("Ret: %f, Offset: %f, Amp: %f, Angle offset: %f\n",ret,fcos->GetParameter(0), fcos->GetParameter(1), fcos->GetParameter(2));
    mg->Add(gr[i]);
  }
  mg->Draw("alp");
  TCanvas *cfit = new TCanvas("cfit","cfit",0,0,1600,900);
  cfit->Divide(2,2);
  TGraph *grf[4];
  //gStyle->SetOptFit(1111);
  TF1 *f0 = new TF1("f0","[0]+[1]*cos(x-[2])",0,7.8);
  TF1 *f1 = new TF1("f1","[0]*abs(sin(x-[1]))",0,7.8);
  for(int i=0;i<3;++i){
    cfit->cd(i+1);
    grf[i] = new TGraph(NN,retard,fit[i]);
    grf[i]->SetMarkerStyle(22);
    grf[i]->SetMarkerSize(0.3);
    grf[i]->Draw("ap");
    if(i==0){
      grf[0]->Fit(f0);
    }
    if(i==1){
      grf[1]->Fit(f1);
    }if(i==2){
      grf[2]->SetTitle("Angle Offset from Fit vs #Delta");
    }
  }
  cfit->cd(4);
  grf[3] = new TGraph(NN,retard,pol);
  grf[3]->SetTitle("Polarization vs #Delta");
  grf[3]->SetMarkerStyle(22);
  grf[3]->SetMarkerSize(0.2);
  grf[3]->Draw("ap");
  grf[3]->Fit(f1);
  return normsq(&V_re_init, &V_im_init);
}
