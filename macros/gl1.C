#include <iostream>
#include "TGraph.h"
#include "TF1.h"

double glcn(double current, int magnet){
  double gl1, gl2, gln=0, fld;
  double cn = current/300.0;
  if(magnet == 1)
    {
      // Moller Quad Q1/MQO1H01/LARGE/new/white
      gl1 = (.0110605+5.33237*cn-.0142794*pow(cn,2)+.259313*pow(cn,3));
      gln = (gl1+0.0058174*pow(cn,4)-0.831887*pow(cn,5));
      fld = gln*10.*5.08/36.5723;
    }
  else if(magnet == 2)
    {
      // Moller Quad Q2/PATSY/MQM1H02/SMALL/RED
      gl1=(0.0196438+5.35443*cn+0.0297273*pow(cn,2)+0.103505*pow(cn,3));
      gln=(gl1-0.0449275*pow(cn,4)-0.211868*pow(cn,5));
      fld=gln*10.*5.08/44.76;
    }
  else if(magnet == 3)
    {
      // Moller Quad Q3/TESSA/MQO1H03/LARGE/BLUE
      gl1=(0.000632446+5.15178*cn-0.00262778*pow(cn,2));
      gl2=(-0.107635*pow(cn,3)+0.00209902*pow(cn,4));
      gln=(gl1+gl2-0.640635*pow(cn,5));
      fld=gln*10.*5.08/36.74 ;
    }
  else if(magnet == 4)
    {
      // Moller Quad Q4/FELICIA/MQO1H03A/LARGE/BLUE
      gl1=(0.0001732+5.2119*cn-0.000732518*pow(cn,2));
      gl2=(-0.133423*pow(cn,3)+0.000618402*pow(cn,4));
      gln=(gl1+gl2-0.647082*pow(cn,5));
      fld=gln*10.*5.08/36.50;
    }
  else if(magnet == 5)
    {
      // Moller Dipole LILLY/MMA1H01/Blue
      fld=(-0.39026E-04+0.027051*current-0.17799E-08*pow(current,2));
    }
  else
    {
      //wrong magnet setup
      fld = 0.0;
    }

  return gln*1e4;
  //return fld;
}
int gl(int mag=1){
  TGraph *gr = new TGraph();
  for(int i=0;i<301;i+=20){
    gr->SetPoint(i,double(i),glcn(double(i),mag));
  }
  gr->SetTitle(Form("Quad %i Pole Field vs. Current", mag));
  gr->SetLineColor(kBlue);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(8);
  gr->Draw("ap");
  return 0;
}
