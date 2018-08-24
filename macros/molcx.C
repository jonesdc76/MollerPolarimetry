#include <iostream>
#include "TMath.h"
#include "TF1.h"

const double me = 0.0005109989461;//mass of electron in GeV
const double alpha = 0.0072973525664;//fine structure constant
//"Calculated cross sections are often given in terms of gigaelectronvolts (GeV),
//via the conversion ħ2c2/GeV2 = 0.3894 mb = 38 940 am2.
//In natural units (where ħ = c = 1), this simplifies to GeV−2 = 0.3894 mb = 38 940 am2."
//see https://en.wikipedia.org/wiki/Barn_(unit)
const double GeV2mb = 0.3894;

//convert beam energy to center of mass energy
//T=kinetic energy of beam in GeV
double Ecom(double T){
  double E = T + me;//convert to total energy
  double Ecm = 2*sqrt(me*(E + me)/2.0);//COM relativistic energy
  return Ecm;//center of mass energy in GeV
}

//Moller scattering cross section in millibarns
//theta1 = lower bound of COM angles in degrees eg. 70
//theta2 = upper bound of COM angles in degrees eg 110
//dphi = range of angles included in phi generation eg. 40
//T = kinetic energy of electron beam in GeV
//see equations 97 to 100 of ./MollerScattering.pdf
double molcx(double theta1=75, double theta2=105, double dphi=40, double T=11){
  double low = cos(theta2*TMath::Pi()/180.0);
  double high = cos(theta1*TMath::Pi()/180.0);
  double Ecm = Ecom(T);//COM relativistic energy
  double coeff = 2*TMath::Pi()*GeV2mb*pow(alpha/Ecm,2)*dphi/180.0;
  coeff *= 0.5; //account for identical particles
  cout<<"Center of mass frame energy: "<<Ecm<<" GeV"<<endl;
  cout<<"Coefficient: "<<coeff<<endl;
  TF1 *f = new TF1("f","pow((x*x+3)/(x*x-1),2)", low, high);
  double cx = coeff;
  cx *= f->Integral(low,high);
  printf("Moller scattering cross section for these parameters: %e mb\n", cx);
  return cx;
}

//Same as molcx but using Bjorken and Drell's Relativistic Quantum Mechanics
double molcxBD(double theta1=75, double theta2=105, double dphi=40, double T=11){
  double low = theta1*TMath::Pi()/180.0;
  double high = theta2*TMath::Pi()/180.0;
  double Ecm = Ecom(T);//COM relativistic energy
  Ecm *= 0.5;//Bjorken and Drell Ecm is per electron
  double coeff = pow(alpha/Ecm,2)/8.0*GeV2mb*2*TMath::Pi()*dphi/180.0;
  coeff *= 0.5; //account for identical particles
  cout<<"Center of mass frame energy: "<<Ecm<<" GeV"<<endl;
  cout<<"Coefficient: "<<coeff<<endl;
  TString fn = Form("sin(x)*((1+pow(cos(x/2.0),4))/pow(sin(x/2.0),4) + 2/pow(sin(x/2)*cos(x/2.0),2)+(1+pow(sin(x/2.0),4))/pow(cos(x/2.0),4))");
  cout<<fn.Data()<<endl;
  TF1 *f = new TF1("f",fn.Data(), low, high);
  double cx = coeff;
  cx *= f->Integral(low, high);
  printf("Moller scattering cross section for these parameters: %e mb\n", cx);
  return cx;
}

//Moller scattering rate in events/sec/microamp for Fe foil
//theta1 = lower bound of COM angles in degrees eg. 70
//theta2 = upper bound of COM angles in degrees eg 110
//dphi = range of angles included in phi generation eg. 40
//T = kinetic energy of electron beam in GeV
//th = thickness of Fe foil in microns
double molrate(double theta1=75, double theta2=105, double dphi=40, double T=11, double th = 1){
  const double rhoFe = 7874;//density of Fe in kg/m^3
  const double avogadro = 6.0221409e26;//atoms per kilomole
  const double atmFe = 55.845;//mass of 1 atom of Fe in amu
  const double qe = 1.60217662e-13;//microcoulomb per electron
  const double Z_Fe = 26;//electrons per atom
  double cx = molcx(theta1, theta2, dphi, T);
  double mrate = 1e-31 * cx * th * 1e-6 / qe * rhoFe / atmFe * avogadro * Z_Fe;
  printf("Moller scattering rate per microamp for %0.2f micron Fe foil: %i events per second\n", th, int(mrate));
  return mrate;
}

//Moller scattering rate in events/sec/microamp for 1 m LH2 target
//theta1 = lower bound of COM angles in degrees eg. 70
//theta2 = upper bound of COM angles in degrees eg 110
//dphi = range of angles included in phi generation eg. 40
//T = kinetic energy of electron beam in GeV
double molrateH(double theta1=75, double theta2=105, double dphi=40, double T=11){
  const double rhoH = 70.8;//density of Fe in kg/m^3
  const double avogadro = 6.0221409e26;//atoms per kilomole
  const double atmH = 1.008;//mass of 1 atom of H in amu
  const double qe = 1.60217662e-13;//microcoulomb per electron
  const double Z_H = 1;//electrons per atom
  double cx = molcx(theta1, theta2, dphi, T);
  double mrate = 1e-31 * cx / qe * rhoH / atmH * avogadro * Z_H;
  printf("Moller scattering rate per microamp for 1 m LH2 taret: %i events per second\n", int(mrate));
  return mrate;
}


//Bjorken and Drell's Moller differential cross section calculation 
double diffmolcxBD(double theta=70.9, double T=0.0157){
  double angle = theta * TMath::Pi() / 180.0;
  double Ecm = Ecom(T);
  Ecm *= 0.5;//Bjorken and Drell Ecm is per electron
  cout<<"Center of mass energy: "<<Ecm<<endl;
  double coeff = pow(alpha/Ecm,2)/8.0*GeV2mb;
  TString fn = Form("((1+pow(cos(x/2.0),4))/pow(sin(x/2.0),4) + 2.0/pow(sin(x/2.0)*cos(x/2.0),2)+(1+pow(sin(x/2.0),4))/pow(cos(x/2.0),4))");
  cout<<fn.Data()<<endl;
  TF1 *f = new TF1("f",fn.Data(), 0,1);
  double diffcx = coeff;
  diffcx *= f->Eval(angle)*1e-3;
  printf("Differential Moller scattering cross section for these parameters: %e barns\n", diffcx);
  return diffcx;
}
