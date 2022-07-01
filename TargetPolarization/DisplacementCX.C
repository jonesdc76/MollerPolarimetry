#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "math.h"
#include "stdio.h"
#include "TPad.h"

const double alpha = 1/137.0, Z = 26, hbarc = 0.197327e-12, M = 55.85 * 931.5;
double m = 0.511, T_d = 0.000040;
double E_beam = 11000;
const double PI = 3.1415927;

Long64_t N_d(double x, double T){
  //number of displaced atoms by initial struck atom is given by the
  //kinetic energy of the atom divided by the displacement energy
  //x = cos(theta) with theta the electron lab scattering angle
  //https://inis.iaea.org/collection/NCLCollectionStore/_Public/46/066/46066650.pdf
  //http://www.nickdelves.co.uk/norgett/norgett/1975_04_23_Calculating_Displacement_Dose_Rates.pdf
  double ke = (T*T/M*(1-x))/(1+T/M*(1-x));
  return Long64_t(ke/(ke<2.5*T_d ? T_d : 2.5*T_d));
}

//Trying to estimate the fractional number of Fe atoms in our target foil displaced in the lattice
//due to beam radiation.
//I attempt using the Mott cross section to get the energy distribution and probability of nuclear recoil
//which I multiply by the number of other atoms displaced by the original recoil N_d() to get the effective
//displacement cross section. I used 40eV for the threshold for displacement using the example of
//http://www.nickdelves.co.uk/norgett/norgett/1975_04_23_Calculating_Displacement_Dose_Rates.pdf
//and integrate from infinity down to threshold in recoil energy or from x=cos(theta) from
//-1 to 0.9999999828 in electron kinematics.
//To convert from nuclear recoil kinetic energy to x use Eq2a from
//https://cds.cern.ch/record/275318/files/CERN-69-22.pdf
double DisplacementCX(double T){
  const int N = 1980;
  double x[N], y[N];
  TF1 *fRuth = new TF1("fRuth", Form("%e*pow(%e/(1-x),2)", PI/2., Z*alpha*hbarc/T), -1, 1);
  TF1 *fMott = new TF1("fMott", Form("%e*pow(%e/(1-x),2)*0.5*(1+x)/(1+(1-x)*%e)",PI/2., Z*alpha*hbarc/T, T/M), -1, 1);

  fMott->SetLineColor(kBlue);
  fMott->Draw();
  fMott->GetXaxis()->SetTitle("cos(#theta)");
  fMott->GetYaxis()->SetTitle("d#sigma/dcos#theta (m^{2})");
  for(int i=0;i<N;++i){x[i] = i*0.001-0.99;y[i]=fRuth->Eval(x[i]);}
  TGraph *gr = new TGraph(N,x,y);
  gPad->SetLogy();
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->Draw("samel");
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  c->SetLogy();
  TF1 *fEp = new TF1("fEp",Form("(%e*(1-x))/(1+%e*(1-x))",T*T/M,T/M),-1,1);
  fEp->GetXaxis()->SetTitle("cos(#theta)");
  fEp->GetYaxis()->SetTitle("Fe Nuclear Recoil Kinetic Energy (MeV)");
  fEp->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",0,0,700,500);
  c2->SetLogy();
  TF1 *fNd = new TF1("fNd",Form("(%e*(1-x))/(1+%e*(1-x))/%e",T*T/M,T/M,T_d),-1,1);
  fNd->GetXaxis()->SetTitle("cos(#theta)");
  fNd->GetYaxis()->SetTitle("# of Displaced Atoms per Recoil Nucleus");
  fNd->Draw();
  TCanvas *c3 = new TCanvas("c3","c3",0,0,700,500);
  c3->SetLogy();
  
  const int n=200000;
  double x1[n], y1[n];
  for(int i=0;i<n;++i){
    x1[i] = i*2.0/double(n)-1.0;
    y1[i] = fMott->Eval(x1[i]);
    //y1[i]= (double)N_d(x[i],T);
    if(i%10==0)printf("%i %e    %e\n", i, x1[i], y1[i]);
  }
  TGraph *gr1 = new TGraph(n,x1,y1);
  //  gPad->SetLogy();
  gr1->SetLineColor(kRed);
  gr1->SetLineWidth(2);
  gr1->Draw("lp");
  TF1 *f_effCX = new TF1("f_effCX",Form("(%e*pow(%e/(1-x),2)" //Rutherford
					"* 0.5*(1+x)/(1+(1-x)*%e))" //Adding in to get to Mott
					"* (%e*(1-x))/(1+%e*(1-x))/%e", //Weight by number of displacements per recoil
					PI/2., Z*alpha*hbarc/T,
					T/M,
					T*T/M, T/M, T_d),
			 -1,1);
  f_effCX->Draw();
  f_effCX->GetXaxis()->SetTitle("cos(#theta)");
  f_effCX->GetYaxis()->SetTitle("Effective Atomic Displacement Cross Section (m^{2})");
  double integ = 0, upper = 0.9999999828, interval = 1e-10, lower = upper - interval;
  int i=0;
  while (lower>-1){
    integ += f_effCX->Integral(lower, upper);
    printf("%i Lower: %e  Upper: %e  Total: %e\n",i, lower, upper, integ);
    ++i;
    interval *=1.2;
    upper = lower;
    lower = upper-interval;
  }
  printf("Mott: %e\n", fMott->Integral(0.99998134224,0.999999253687));
  return f_effCX->Integral(-1,0.9999999828);
}
