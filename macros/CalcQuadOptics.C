//*****************************
//Author:  Donald Jones      *
//Contact: jonesdc@jlab.org  *
//Date:    02/10/2023        *
//*****************************
//*************************************************************************************************************
// Program to calculate optics tracks for Moller scattering through the 4 quadrupole magnets and the dipole  *
// aperture in the Hall A Moller polarimeter spectrometer. If an ideal track is determined, the program will *
// perform a Chi-square fit to determine the quadrupole currents that  minimize the deviations from the      *
// ideal track.                                                                                              *
// CalcQuadOptics() is the main program and it calls the function getX() where the optics calculations are   *
// performed.                                                                                                *
// Expects data file quad_currents.dat with quadrupole currents.                                             *
//*************************************************************************************************************

#include "TF1.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TLine.h"

const double PI = 3.1415927;

//TGTZ is a parameter for moving the target z position. Set this to 0 for the setup prior to 2023
//and to 0.3 for 2024 and after when the target was moved upstream 30 cm to accommodate higher energy.
const double TGT_Z = 0.30;//

//Quadrupole and dipole z-position information relative to target
//center "_Z", effective length "_L", z of upstream edge "_u" and z of downstream edge "_d"
const double Q1_Z = 0.7519 + TGT_Z, Q1_L = 0.3658, Q1_u = Q1_Z - Q1_L/2.0, Q1_d = Q1_Z + Q1_L/2.0;
const double Q2_Z = 1.4046 + TGT_Z, Q2_L = 0.4477, Q2_u = Q2_Z - Q2_L/2.0, Q2_d = Q2_Z + Q2_L/2.0;
const double Q3_Z = 2.0908 + TGT_Z, Q3_L = 0.3674, Q3_u = Q3_Z - Q3_L/2.0, Q3_d = Q3_Z + Q3_L/2.0;
const double Q4_Z = 2.7459 + TGT_Z, Q4_L = 0.3650, Q4_u = Q4_Z - Q4_L/2.0, Q4_d = Q4_Z + Q4_L/2.0;
const double Dip_Z = 4.234 + TGT_Z, Dip_L = 1.645, Dip_u = Dip_Z- Dip_L/2.0, Dip_d = Dip_Z+ Dip_L/2.0;


const double SLIT_X = 0.045;//Horizontal distance from beam to dipole slit center
const double me = 0.000511;//electron mass in GeV/c^2

//A few bad practice global defaults
int Npasses = 4;//default number of beam passes through the accelerator
double E_beam = 2.08183*double(Npasses)+0.1185; //default electron beam energy in GeV
double p_beam = sqrt(E_beam * E_beam - me*me);//default beam momentum in GeV/c
bool verbose = false;

using namespace std;

//*****************************************************************************************************
// getX() calculates the transverse horizontal distance (X) from the unscattered beam trajectory     *
// for Moller scattered electrons through the four quadrupoles in the Hall A Moller polarimeter.     *
// The format of this function is set by the requirements of the TF1 fitting function that uses      *
// it to find optics solutions.                                                                      *
//                                                                                                   *
// The quadrupole field parameterization as a function of current comes from the following technote: *
// https://github.com/jonesdc76/MollerPolarimetry/blob/master/quads/QuadrupoleInfo.pdf               *
//*****************************************************************************************************
//												     *
// Arguments: *x   - array of the dimension variables. Since this is a function of Z distance along  *
//                   the beam only x[0] is used.						     *
//            *par - array of parameters that can be set or varied in the TF1 fit function.          *
//                   Parameters 0-3 are the four quadrupole currents which are input as fractions of *
//                   the maximum current of 300A.                                                    *
//                   Parameter 4 (fifth parameter) is reserved for the COM scattering angle in       *
//                   degrees and is not allowed to vary in the fit.                                  *
//                                                                                                   *
// Return:    horizontal transverse distance in meters of the scattered electron                     *
//                                                                                                   *
// Note that for optics with beam energies above 5 pass requires changing the quad strengths and/or  *
// target Z position. This is currently implemented by doubling the quad strengths above 5 pass      *
// using the variable "optics22GeV". Effective doubling of the quad strengths is possible per Jay    *
// with a combination of a new design with stronger field and greater length. Changing the constant  *
// TGTZ to a larger value will move the target further upstream.                                     * 
//*****************************************************************************************************

double getX(double *x, double *par){
  double z_pos = x[0];//single independent variable dimension z
  if(z_pos < 0 ){
    printf("Invalid z-position %f. Must be positive.\n", z_pos);
    return 0;
  }
  
  //Define some variables for both COM and lab frames.
  ////////////////////////////////////////////////////
  double scatt_angle_cm = par[4] / 180.0 * PI;//5th paramter gives COM scattering angle in degrees.
  double E_cm = sqrt(me*(E_beam + 2*me)/2.0);//CM energy of single electron after scattering in GeV
  double p_cm = sqrt(E_cm*E_cm-me*me);//CM momentum in GeV/c
  double gamma = sqrt(p_beam/me/2.0);//high energy approximation gamma of relative velocity between lab and CM
  double beta = sqrt(1-pow(0.5*E_beam/me,-2));//lab frame scattered electron beta v/c
  double vel = p_beam/(E_beam+me);//relative velocity between COM and lab frames in units of c
  double tan_theta_lab = sin(scatt_angle_cm)/(gamma*(cos(scatt_angle_cm)+vel*E_cm/p_cm));
  //tan_theta_lab = p_cm/(gamma*vel*E_cm);

  //Initial vector at scattering vertex [x=0 , tan(theta_lab)]
  ////////////////////////////////////////////////////////////
  TVectorD V_x = TVectorD(2);
  V_x[0] = 0;
  V_x[1] = tan_theta_lab;
  double z = z_pos <= Q1_u ? z_pos : Q1_u;

  
  //Drift between target and quad 1
  /////////////////////////////////////////////
  TMatrixD M_Drift(2,2);
  M_Drift(0,0) = 1;
  M_Drift(0,1) = z;
  M_Drift(1,0) = 0;
  M_Drift(1,1) = 1;

  //Vector at entrance to quad 1 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Drift;
  if(verbose) printf("Q1u: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
    
  //(De)Focussing in quad 1
  /////////////////////////////////////////////
  z = z_pos <= Q1_d ? z_pos : Q1_d;
  double c[6] = {0.007993, 5.233, 0.02832, 0.2878, -0.03376, -0.8568};
  double k = c[0] + c[1]*par[0]+ c[2]*pow(par[0],2) + c[3]*pow(par[0],3) + c[4]*pow(par[0],4) + c[5]*pow(par[0],5);
  k *= 0.2998 / (E_beam/2. * beta * Q1_L);

  TMatrixD M_Q(2,2);
  double d = z - Q1_u;
  double sqrt_k = sqrt(abs(k));
  M_Q(0,0) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);
  M_Q(0,1) = k < 0 ? 1/sqrt_k*sinh(sqrt_k * d) : 1/sqrt_k*sin(sqrt_k * d);
  M_Q(1,0) = k < 0 ? sqrt_k*sinh(sqrt_k * d) : -sqrt_k*sin(sqrt_k * d);
  M_Q(1,1) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);

  //Vector at exit to quad 1 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Q;
  if(verbose) printf("Q1d: position %f   angle %f\n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
  
  //Drift between quad 1 and quad 2
  //////////////////////////////////////////////
  z = z_pos <= Q2_u ? z_pos : Q2_u;
  M_Drift(0,0) = 1;
  M_Drift(0,1) = z - Q1_d;
  M_Drift(1,0) = 0;
  M_Drift(1,1) = 1;

  //Vector at entrance to quad 2 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Drift;
     if(verbose) printf("Q2u: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
    
  //(De)Focussing in quad 2
  //////////////////////////////////////////////
  z = z_pos <= Q2_d ? z_pos : Q2_d;
  d = z - Q2_u;
  c[0] = 0.0234; c[1] = 5.3544; c[2] = 0.0135, c[3] = 0.1038, c[4] = -0.0318, c[5] = -0.2121;
  k = c[0] + c[1]*par[1]+ c[2]*pow(par[1],2) + c[3]*pow(par[1],3) + c[4]*pow(par[1],4) + c[5]*pow(par[1],5);
  k *= 0.2998 / (E_beam/2. * beta * Q2_L);

  sqrt_k = sqrt(abs(k));
  M_Q(0,0) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);
  M_Q(0,1) = k < 0 ? 1/sqrt_k*sinh(sqrt_k * d) : 1/sqrt_k*sin(sqrt_k * d);
  M_Q(1,0) = k < 0 ? sqrt_k*sinh(sqrt_k * d) : -sqrt_k*sin(sqrt_k * d);
  M_Q(1,1) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);

  //Vector at exit to quad 2 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Q;
  if(verbose) printf("Q2d: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
  
  //Drift between quad 2 and quad 3
  //////////////////////////////////////////////
  z = z_pos <= Q3_u ? z_pos : Q3_u;
  M_Drift(0,0) = 1;
  M_Drift(0,1) = z - Q2_d;
  M_Drift(1,0) = 0;
  M_Drift(1,1) = 1;

  //Vector at entrance to quad 3 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Drift;
  if(verbose) printf("Q3u: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
    
  //(De)Focussing in quad 3
  //////////////////////////////////////////////
  z = z_pos <= Q3_d ? z_pos : Q3_d;
  d = z - Q3_u;
  c[0] = -0.02298; c[1] = 5.3; c[2] = 0.1018, c[3] = 0.163, c[4] = -0.05161, c[5] = -0.2734;
  k = c[0] + c[1]*par[2]+ c[2]*pow(par[2],2) + c[3]*pow(par[2],3) + c[4]*pow(par[2],4) + c[5]*pow(par[2],5);
  k *= 0.2998 / (E_beam/2. * beta * Q3_L);

  sqrt_k = sqrt(abs(k));
  M_Q(0,0) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);
  M_Q(0,1) = k < 0 ? 1/sqrt_k*sinh(sqrt_k * d) : 1/sqrt_k*sin(sqrt_k * d);
  M_Q(1,0) = k < 0 ? sqrt_k*sinh(sqrt_k * d) : -sqrt_k*sin(sqrt_k * d);
  M_Q(1,1) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);

  //Vector at exit to quad 3 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Q;
  if(verbose) printf("Q3d: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
  
  //Drift between quad 3 and quad 4
  //////////////////////////////////////////////
  z = z_pos <= Q4_u ? z_pos : Q4_u;
  M_Drift(0,0) = 1;
  M_Drift(0,1) = z - Q3_d;
  M_Drift(1,0) = 0;
  M_Drift(1,1) = 1;

  //Vector at entrance to quad 4 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Drift;
  if(verbose) printf("Q4u: position %f   angle %f \n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];

      
  //(De)Focussing in quad 4
  //////////////////////////////////////////////
  z = z_pos <= Q4_d ? z_pos : Q4_d;
  d = z - Q4_u;
  c[0] = -0.01803; c[1] = 5.07; c[2] = 0.08835, c[3] = 0.02901, c[4] = -0.04965, c[5] = -0.7211;
  k = c[0] + c[1]*par[3]+ c[2]*pow(par[3],2) + c[3]*pow(par[3],3) + c[4]*pow(par[3],4) + c[5]*pow(par[3],5);
  k *= 0.2998 / (E_beam/2. * beta * Q4_L);
 		
  sqrt_k = sqrt(abs(k));
  M_Q(0,0) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);
  M_Q(0,1) = k < 0 ? 1/sqrt_k*sinh(sqrt_k * d) : 1/sqrt_k*sin(sqrt_k * d);
  M_Q(1,0) = k < 0 ? sqrt_k*sinh(sqrt_k * d) : -sqrt_k*sin(sqrt_k * d);
  M_Q(1,1) = k < 0 ? cosh(sqrt_k * d) : cos(sqrt_k * d);

  //Vector at exit to quad 4 [x, tan(theta)]
  //////////////////////////////////////////////
  V_x *= M_Q;
  if(verbose) printf("Q4d: position %f   angle %f \n\n", V_x[0], V_x[1]);
  if(z_pos == z)return V_x[0];
  
  //Drift after quad 4
  //////////////////////////////////////////////
  z = z_pos;
  M_Drift(0,0) = 1;
  M_Drift(0,1) = z - Q4_d;
  M_Drift(1,0) = 0;
  M_Drift(1,1) = 1;

  //Vector at downstream end of dipole [x, tan(theta)]
  ////////////////////////////////////////////////////
  V_x *= M_Drift;
  return V_x[0];
  
}



//******************************************************************************************************************
// Uses electron trajectories from optics calculated in getX() to find paths through the Moller spectrometer.     *
// Only calculates horizontal distance X from beam center to determine optimal path through the symmetric         *
// slits in the dipole.                                                                                           *
// Reads initial quadrupole current setpoints from data file "quad_currents.dat" and can simply plot trajectories *
// for these optics or if an ideal trajectory is defined, it can minimize the deviation from ideal with a Chi-    *
// squared minimization to get the current setpoints that produce the closest to the ideal trajectory.            *
//******************************************************************************************************************
//                                                                                                                *
// Arguments:  Npass       - number of passes for electron beam through accelerator gives electron beam energy    *
//             deg_range   - number of degrees you want to include around 90 degrees COM scattering. For example  *
//                           to calculate ray traces for 90+/-10 choose deg_range = 10                            *
//             optimize    - choose true if you want to fit to find the optimal currents to obtain as close to    *
//                           the ideal electron trajectory as possible.                                           *
//******************************************************************************************************************


int CalcQuadOptics(int Npass = 1, double deg_range = 10, bool optimize = false){

  verbose = false;
  
  //Set up plot style
  /////////////////////////////////////////////
  gStyle->SetStatX(0.35);
  gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.125);
  gStyle->SetPadRightMargin(0.05);
  
  //Set up drawing canvas
  /////////////////////////////////////////////
  TCanvas *c1 = new TCanvas("c1","Hall A Moller Quad Optics",1200,800);
  c1->SetTickx();
  c1->SetTicky();
  
  //Allow from 1 to 5 passes (up to 11 GeV)
  /////////////////////////////////////////////
  if(Npass > 0 && Npass < 6){
    Npasses = Npass;
  }else{
    printf("\n!!The optics of the Moller spectromter for the 11 GeV era do not allow for higher energies."
	   " Calculating 5 pass optics instead!!\n\n");
    Npasses = 5;
  }

  //Set up optical parameters to solve for
  /////////////////////////////////////////////
  E_beam = 2.08183*double(Npasses)+0.1185; //Electron beam energy in GeV Spring 2023
  p_beam = sqrt(E_beam * E_beam - me*me);  //beam momentum in GeV/c


  //Grab optics values from file
  ///////////////////////////////
  ifstream datafile;
  datafile.open("quad_currents.dat");
  char n[100], q1[100], q2[100], q3[100], q4[100];
  double q_curr[4];int np;
  while ( !datafile.eof() ){
    datafile >>n>>q1>>q2>>q3>>q4;
    if(atoi(n)==Npasses){
      np = atoi(n);
      q_curr[0] = atof(q1);
      q_curr[1] = atof(q2);
      q_curr[2] = atof(q3);
      q_curr[3] = atof(q4);
      printf("Quad currents at %d-pass:  %f %f %f %f\n",
	     np, q_curr[0], q_curr[1], q_curr[2], q_curr[3]);
      break; 
    }
  }
  datafile.close();


  //Main plot with all the traces and optics information
  ///////////////////////////////////////////////////////
  const int N = 100;
  double ytop = 0.06;//top of plot Y axis
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(Form("Optics of 90#circ#pm%i#circ CM Moller Scattered Electrons in Hall A Moller Polarimeter",
		    int(deg_range)));


  //Set up areas to shade on plot indicating where the quadrupoles and dipole are.
  /////////////////////////////////////////////////////////////////////////////////
  double dist_Z = Dip_d;//terminal of optical trace (X-axis upper limit) usually set to downstream end of dipole
  double xin = SLIT_X-0.01, xout = SLIT_X+0.01;//distance from beam center to inner and outer walls of dipole slit
  double dipz[5]={Dip_u,Dip_d,Dip_d,Dip_u,Dip_u}, dipx1[5]={0,0,xin,xin,0}, dipx2[5]={xout,xout,ytop,ytop,xout};
  TGraph *grdip1 = new TGraph(5, dipz, dipx1);
  TGraph *grdip2 = new TGraph(5, dipz, dipx2);
  grdip1->SetFillColor(1);
  grdip1->SetFillStyle(1001);
  grdip2->SetFillColor(1);
  grdip2->SetFillStyle(1001);
  double quadu[4] = {Q1_u,Q2_u,Q3_u,Q4_u}, quadd[4] = {Q1_d,Q2_d,Q3_d,Q4_d};
  TGraph *grq[4];
  for(int i=0;i<4;++i){
    double ztmp[5] = {quadu[i],quadd[i],quadd[i],quadu[i],quadu[i]};
    double xtmp[5] = {0,0,ytop,ytop,0};
    grq[i] = new TGraph(5,ztmp,xtmp);
    grq[i]->SetFillColor(38);
    grq[i]->SetFillStyle(1001);
    grq[i]->GetXaxis()->SetLimits(0,dist_Z);
    grq[i]->GetYaxis()->SetLimits(0,ytop);
    mg->Add(grq[i],"F");
  }

  //Draw in edges of beam pipe
  /////////////////////////////////////////////
  double zpipe[5] = {0,Q4_d+0.1,Q4_d+0.1,0,0}, xpipe[5]={0.05,0.05,0.0505,0.0505,0.05};
  TGraph *grpipe1 = new TGraph(5,zpipe,xpipe);
  grpipe1->SetFillColor(kBlack);
  grpipe1->SetFillStyle(1011);
  mg->Add(grpipe1,"F");

  zpipe[0] = Q4_d+0.1; zpipe[1]= Q4_d+0.1; zpipe[2] = Q4_d+0.08; zpipe[3] = Q4_d+0.08; zpipe[4] = Q4_d+0.1;
  xpipe[0] = 0.05; xpipe[1] = xout; xpipe[2]=xout; xpipe[3]=0.05; xpipe[4] = 0.05;
  TGraph *grpipe2 = new TGraph(5,zpipe,xpipe);
  grpipe2->SetFillColor(kBlack);
  grpipe2->SetFillStyle(1011);
  mg->Add(grpipe2,"F");

  zpipe[0] = Q4_d+0.08; zpipe[1]= Dip_u; zpipe[2] = Dip_u; zpipe[3] = Q4_d+0.08; zpipe[4] = Q4_d+0.08;
  xpipe[0] = xout; xpipe[1] = xout; xpipe[2]=xout+0.0005; xpipe[3]=xout+0.0005; xpipe[4] = xout;
  TGraph *grpipe3 = new TGraph(5,zpipe,xpipe);
  grpipe3->SetFillColor(kBlack);
  grpipe3->SetFillStyle(1011);
  mg->Add(grpipe3,"F");

  //If desired, set ideal electron trajectory
  /////////////////////////////////////////////
  double x[N+1], xe[N+1], z[N+1], ze[N+1];
  double slope = 0.0020, offset = (xin+0.006);
  //Straight sloped line from target to Q4 with slope change after Q4 through dipole slit
  for(int i=0; i<N+1; ++i){
    z[i] = i * dist_Z/double(N);
    ze[i] = 0;
    x[i] = (z[i] > Q4_Z ? offset + slope*(z[i]-Q4_Z) : offset * z[i]/Q4_Z);
    xe[i] = z[i]<1? 0.001 : 0.001/z[i];
  }
  TGraphErrors *grIdeal = new TGraphErrors(N+1, z, x, ze, xe);
  grIdeal->SetMarkerStyle(8);
  grIdeal->SetMarkerSize(0.7);
  grIdeal->SetMarkerColor(kGreen+2);
  grIdeal->SetLineColor(kGray);
  
  //Five parameter function. First 4 parameters (quad currents) are fit parameters.
  //Quad currents are input as a fraction of maximum current which is 300 A.
  //5th parameter is fixed and is the COM scattering angle in degrees.
  /////////////////////////////////////////////////////////////////////////////////
  TF1 *f = new TF1("f", getX, 0, dist_Z, 5);
  f->SetLineColor(kBlack);
  f->SetLineStyle(2);
  f->SetLineWidth(5);
  f->SetParNames("Q1 Current", "Q2 Current", "Q3 Current", "Q4 Current", "#theta_{CM}");
  
  //Choose quadrupole currents either as initialization values for the fit or from
  //prior knowledge. Currents are fractions of maximum current which is 300 A.
  /////////////////////////////////////////////////////////////////////////////////
  f->SetParameters(q_curr[0]/300., q_curr[1]/300., q_curr[2]/300., q_curr[3]/300., 90);
  /* if(Npasses==1) f->SetParameters( 0.333, 0.0786, 0.000, 0.138, 90);//1-pass optics */
  /* if(Npasses==2) f->SetParameters(-0.137, 0.324,  0.140, 0.138, 90);//2-pass optics */
  /* if(Npasses==3) f->SetParameters(-0.611, 0.202,  0.281, 0.311, 90);//3-pass optics */
  /* if(Npasses==4) f->SetParameters(-0.611,-0.746,  0.620, 0.625, 90);//4-pass optics */
  /* if(Npasses==5) f->SetParameters(-0.682,-0.823,  0.498, 0.724, 90);//5-pass optics */
  f->FixParameter(4, 90);//Don't use 5th parameter in fit
  for(int i=0;i<4;++i)f->SetParLimits(i,-1,1); //Limit currents to physical limits.

  //Either fit to find currents giving closest to  ideal trajectory or plot
  //trajectory of chosen optics solution.
  //////////////////////////////////////////////////////////////////////////////////
  if(!optimize){
    for(int i=0;i<4;++i) f->FixParameter(i,f->GetParameter(i));
  }else{
    mg->Add(grIdeal,"lp");
  }
  grIdeal->Fit(f,"r");
  mg->Draw("A");

  //Now draw trajectories in range around theta_cm=90 degrees
  //////////////////////////////////////////////////////////////////////////////////
  f->FixParameter(4, 90-deg_range);
  for(int i=0;i<N;++i){
    z[i] = i*Dip_d/double(N-1);
    x[i] = f->Eval(z[i]);
  }
  TGraph *grLower = new TGraph(N,z,x);
  grLower->SetLineColor(12);
  grLower->SetLineWidth(2);
  mg->Add(grLower,"c");
  f->FixParameter(4, 90);
  for(int i=0;i<N;++i){
    z[i] = i*Dip_d/double(N-1);
    x[i] = f->Eval(z[i]);
  }
  TGraph *gr90 = new TGraph(N,z,x);
  gr90->SetLineColor(1);
  gr90->SetLineWidth(5);
  gr90->SetLineStyle(10);
  mg->Add(gr90,"c");
  f->FixParameter(4, 90+deg_range);
  for(int i=0;i<N;++i){
    z[i] = i*Dip_d/double(N-1);
    x[i] = f->Eval(z[i]);
  }
  TGraph *grUpper = new TGraph(N,z,x);
  grUpper->SetLineColor(12);
  grUpper->SetLineWidth(2);
  mg->Add(grUpper,"c");
  mg->Add(grdip1,"F");
  mg->Add(grdip2,"F");


  //Label and clean up the plot
  //////////////////////////////////////////////////////////////////////////////////
  mg->GetXaxis()->SetTitle("Z Distance along Beam from Target (m)");
  mg->GetYaxis()->SetTitle("X Horizontal Distance from Beam (m)");
  mg->GetXaxis()->SetLimits(0,dist_Z);
  mg->GetXaxis()->SetRangeUser(0,dist_Z);
  mg->GetYaxis()->SetRangeUser(0,ytop);
  gPad->RedrawAxis();
  gPad->RedrawAxis("G");
  gPad->Update();


  //Print quadrupole current settings on plot
  //////////////////////////////////////////////////////////////////////////////////  
  TPaveText *pt1 = new TPaveText(Q1_Z-0.2,0.85*ytop,Q1_Z+0.2,0.95*ytop);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->SetFillStyle(0);
  pt1->AddText("Q1");
  pt1->AddText(Form("%+4.1fA",f->GetParameter(0)*300));
  pt1->Draw();
  TPaveText *pt2 = new TPaveText(Q2_Z-0.2,0.85*ytop,Q2_Z+0.2,0.95*ytop);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->SetFillStyle(0);
  pt2->AddText("Q2");
  pt2->AddText(Form("%+4.1fA",f->GetParameter(1)*300));
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(Q3_Z-0.2,0.85*ytop,Q3_Z+0.2,0.95*ytop);
  pt3->SetBorderSize(0);
  pt3->SetFillColor(0);
  pt3->SetFillStyle(0);
  pt3->AddText("Q3");
  pt3->AddText(Form("%+4.1fA",f->GetParameter(2)*300));
  pt3->Draw();
  TPaveText *pt4 = new TPaveText(Q4_Z-0.2,0.85*ytop,Q4_Z+0.2,0.95*ytop);
  pt4->SetBorderSize(0);
  pt4->SetFillColor(0);
  pt4->SetFillStyle(0);
  pt4->AddText("Q4");
  pt4->AddText(Form("%+4.1fA",f->GetParameter(3)*300));
  pt4->Draw();
  TPaveText *pt5 = new TPaveText(Dip_Z-0.4,0.025,Dip_Z+0.4,0.030);
  pt5->SetBorderSize(0);
  pt5->SetFillColor(1);
  pt5->SetTextColor(kWhite);
  pt5->AddText("Dipole");
  pt5->Draw();
  TPaveText *pt6 = new TPaveText(Q1_Z/2.0-0.3,0.049,Q1_Z/2.0+0.3,0.0515);
  pt6->SetBorderSize(0);
  pt6->SetFillColor(0);
  pt6->SetTextColor(1);
  pt6->AddText("Beam Pipe");
  pt6->Draw();
  TPaveText *pt7 = new TPaveText(0.15,0.053,Q1_u-0.04,0.058);
  //pt7->SetBorderSize(0);
  pt7->SetShadowColor(0);
  pt7->SetFillColor(0);
  pt7->SetTextColor(1);
  pt7->AddText(Form("E_{beam}=%0.1f GeV",E_beam));
  pt7->Draw();

  
  return 0;
}
