/*
Don Jones, June 2017
...........................

Calculate magnetization corrections as a function of temperature and B-field
Follows treatment in Equations 4 and 7 in 
"Spin-waves in nickel, iron and yttrium-iron garnet", Pauthenet, March 1982

Note: Start with room-temperature measurements of magnetization at 1 Telsa
      then only need to calculate the small corrections from increases in
      temperature and applied B-field. These corrections are given by dM/dH
      in Eq 7, but since the values of A and B are not given for room 
      temperature we will have to differentiate Eq 9 w.r.t. H and T. The 
      internal field H_i in these equations is the total field B including the
      applied field plus magnetization but above the saturation magnetization
      dM/dB=dM/dH_i=dM/dH. Slopes are assumed to be linear over the small ranges
      of temperature and magnetization such that (delta M) = dM/dH*(delta H) 
      and (delta M) = dM/dT*(delta T)

Arguments: meas_T is the temperature in Kelvin at which the magnetization 
           is known from measurement(assumed to be room temperature)

           meas_Hi is the internal field in Oersteds at which the magnetization 
	   was measured (Above saturation, the relation is a nearly constant 
	   offset of H = H_i+7.177 Oe for Fe. Below saturation the relationship
	   is not so simple and this program is insufficient. If the applied 
	   field is assumed to be 10000 Oersteds = 1 Telsa, then the internal 
	   field H_i = 2823 Oe.)

Return:  total correction to be ADDED to starting magnetization to correct
         for effects from temperature and applied field increases

NOTE on translating from internal field H_i to applied field H.
-------------------------------------------------------------------------------
      The Pauthenet paper calculates corrections in terms of internal field H_i.
      The internal field H_i is the difference between the applied field & the  
      demagnetizing field (https://en.wikipedia.org/wiki/Demagnetizing_field)
      H_i = H - gM (emu system) where g is the demagnetizing factor and M is  
      magnetization. For spheres in a uniform applied field H, g = 4pi/3 giving       H_i = H - 4piM/3 (see eq 5.98 Jackson 2nd edition). Pauthenet used 5 mm  
      diameter spheres in this measurement and since we generally want to talk 
      in terms of applied field we can use H = H_i+4piM/3 to translate his 
      corrections to be in terms of the applied field H. Above saturation, the 
      magnetization is approximately constant. For iron the saturation 
      magnetization is given as 217.6+/-0.1 emu/g  at T=293 K and H=10 kOe in 
      Crangle and Goodman(1970)which is equivalent to 4pirho(217.6)=21.53 kOe 
      (rho is density in g/cm^3) so 4piM/3= 7.177 kOe. Thus, to translate from 
      internal field to H you add 7.177 kOe. My default field is H=10 kOe which
      gives H_i=10.000-7.177 kOe and my default temperature is 293 K.
*/ 
#include <iostream>
#include <cstdio>
const double a3_2 = 307e-6, a5_2 = 22.8e-8, b = 1.378e-4;
const double M43rdsPi = 7177; //magnetization of Iron in Oersted. Add this to
                              // H_i to get applied field H above saturation.
double CorrectionsMagnetization(double dT, double dH, //in Kelvin and Oersted
				double meas_T = 293, double meas_Hi = 2823){

  double dM = 0, M_start, M_final;
  //calculate 307e-6*T^(3/2)*F(3/2, 1.378*Hi/T) from equation (9) at meas_T
  //and meas_Hi
  double s = 3.0/2.0, term1 = 0, term2 = 0, term3 = 0;
  for(int i=1;i<10000;++i){
    term1 += a3_2*pow(meas_T,s)*pow(i,-s)*exp(-i*b*meas_Hi/meas_T);
    term1 -= a3_2*pow(meas_T+dT,s)*pow(i,-s)*exp(-i*b*(meas_Hi+dH)/(meas_T+dT));
  }
  s = 5.0/2.0;
  for(int i=1;i<10000;++i){
    term2 += a5_2*pow(meas_T,s)*pow(i,-s)*exp(-i*b*meas_Hi/meas_T);
    term2 -= a5_2*pow(meas_T+dT,s)*pow(i,-s)*exp(-i*b*(meas_Hi+dH)/(meas_T+dT));
    
  }
  //my linear parameterization of Chi term using fit to values in Table 1
  double intercept = 3.644e-6, slope = 5.0434e-10;
  term3 += intercept*dH + slope * (dT*meas_Hi + (meas_T+dT)*dH);
  cout<<"Term 1:"<<term1<<endl;
  cout<<"Term 2:"<<term2<<endl;
  cout<<"Term 3:"<<term3<<endl;
  dM = term1+term2+term3;
  return dM;
}
