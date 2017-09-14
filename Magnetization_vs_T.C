{
//This graph comes from Crangle and Goodman 1970 "The Magnetization of Pure Iron and Nickel"
  bool useFe = 0;
  const int N=40;
  const double conversion = 0.00999921;//convert from emu/g to emu/atom
  const double CurieTempFe = 1044.0;//Curie temperature of Iron in Kelvin
  const double CurieTempNi = 631.0;//Curie temperature of Nickel in Kelvin
  const double satMag0KFe = 221.71;
  const double satMag0KNi = 58.58;
  double absMagFe[N] = {221.7,221.7,221.6,221.5,221.3,221.1,220.9,220.6,220.3,219.9,219.6,219.1,218.7,218.2,217.6,217.2,216.7,216.1,215.5,214.7,213.8,212.9,212.0,210.9,209.8,208.4,207.1,205.7,204.2,202.6,200.8,199.0,196.9,195.0,192.7,190.3,187.6,184.8,181.5,177.9};
  double absMagNi[N] = {58.58,58.56,58.52,58.48,58.43,58.36,58.3,58.22,58.14,58.06,57.95,57.82,57.68,57.54,57.39,57.22,57.04,56.84,56.61,56.38,56.11,55.81,55.51,55.18,54.68,54.29,53.9,53.44,52.99,52.52,52.03,51.56,51.08,50.57,50.02,49.32,48.53,47.62,46.5,45.17};
  double reducedTemp[N], tempFe[N], reducedMagFe[N], tempNi[N], reducedMagNi[N];
  for(int i=0;i<N;++i){
    reducedTemp[i] = (i*0.02);
    reducedMagFe[i] = absMagFe[i]/satMag0KFe;
    reducedMagNi[i] = absMagNi[i]/satMag0KNi;
    tempFe[i] = reducedTemp[i]*CurieTempFe;
    tempNi[i] = reducedTemp[i]*CurieTempNi;
    if(useFe)
      printf("%0.2f  %0.4f | %0.1f  %0.1f\n",reducedTemp[i], reducedMagFe[i],tempFe[i], absMagFe[i]);
    else
      printf("%0.2f  %0.4f | %0.1f  %0.1f\n",reducedTemp[i], reducedMagNi[i],tempNi[i], absMagNi[i]);
  }
  TF1 *f = new TF1("f","pol1",270,340);
  TGraph *gr;
  if(useFe)
    gr = new TGraph(N, tempFe, absMagFe);
  else
    gr = new TGraph(N, tempNi, absMagNi);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.8);
  gr->SetLineColor(kRed);
  gr->Fit(f,"r");
  TPaveText *pt = new TPaveText(0.5,0.8,0.89,0.89,"ndc");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText("Temperature Sensitivity at 294 K");
  double sens = f->GetParameter(1)/(f->GetParameter(0)+ 294*f->GetParameter(1));
  pt->AddText(Form("%0.4f\%/^{o}C", sens*100));
  if(useFe)
    gr->SetTitle("Temperature Dependence of Saturation Magnetization for Pure Iron");
  else
    gr->SetTitle("Temperature Dependence of Saturation Magnetization for Pure Nickel");
  gr->GetXaxis()->SetTitle("Temperature (K)");
  gr->GetYaxis()->SetTitle("Saturation Magnetization (emu/g)");
  gr->Draw("ap");
  pt->Draw();
  if(0){
    TGraph *grRedFe = new TGraph(N, reducedTemp, reducedMagFe);
    grRedFe->SetMarkerColor(kBlue);
    grRedFe->SetMarkerStyle(8);
    grRedFe->SetMarkerSize(0.8);
    grRedFe->SetLineColor(kRed);
    grRedFe->SetTitle("Temperature Dependence of Saturation Magnetization for Pure Iron(Blue) and Nickel(Red)");
    grRedFe->GetXaxis()->SetTitle("Reduced Temperature");
    grRedFe->GetYaxis()->SetTitle("Reduced Saturation Magnetization");
    grRedFe->Draw("ap");
    TGraph *grRedNi = new TGraph(N, reducedTemp, reducedMagNi);
    grRedNi->SetMarkerColor(kRed);
    grRedNi->SetMarkerStyle(8);
    grRedNi->SetMarkerSize(0.8);
    grRedNi->Draw("samep");
  }
  double T_chi[11] = {4.21,24.79,51.05,75.34,100.51,131.55,165.81,197.45,226.34,254.53,286.41};
  double chi[11] = {3.61,3.71,3.70,3.76,3.67,3.64,3.66,3.66,3.79,3.83,3.82};
  for(int i=0;i<11;++i)chi[i]*=1e-6;
  TGraph *grChi = new TGraph(11,T_chi,chi);
  grChi->SetMarkerColor(kRed);
  grChi->SetMarkerStyle(8);
  //  grChi->Draw("alp");
}

