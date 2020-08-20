#include <iostream>
#include <TGraphErrors.h>
#include <TString.h>
#include <TF1.h>
  const double pi = 3.1415926;
  double theta = pi/3.0;//angle of birefringent optic w.r.t. the x-axis
  double Delta = 0.02;//retardance in radians
  double Re_x1(double theta, double Delta){
    return (pow(cos(theta),2)+pow(sin(theta),2)*cos(Delta)-cos(theta)*sin(theta)*sin(Delta))/sqrt(2);
  }
  double Im_x1(double theta, double Delta){
    return (pow(sin(theta),2)*sin(Delta)+cos(theta)*sin(theta)-cos(theta)*sin(theta)*cos(Delta))/sqrt(2);
  }
  double Re_x2(double theta, double Delta){
    return (cos(theta)*sin(theta)-cos(theta)*sin(theta)*cos(Delta)-pow(cos(theta),2)*sin(Delta))/sqrt(2);
  }
  double Im_x2(double theta, double Delta){
    return (cos(theta)*sin(theta)*sin(Delta)+pow(sin(theta),2)+pow(cos(theta),2)*cos(Delta))/sqrt(2);
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
int birefringence(){
  const int N=18;
  TF1 *fTsq[N];
  for(int i=N-1;i>=0;--i){
    fTsq[i] = new TF1(Form("fTsq%i",N),func(theta, 0.1*i).Data(),0,2*3.1415926);
    fTsq[i]->SetLineColor(kGreen+i);
    fTsq[i]->Draw((i==N-1?"":"same"));
  }
  return 0;
}
