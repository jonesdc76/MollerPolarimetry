{
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetOptFit(1111);
  const int N=8;
  double in[N] = {5.2316,5.2533,5.2157,5.2422,5.2063,5.2786,5.2824,5.2625};
  double inE[N] = {0.0115,0.0116,0.0101,0.0119,0.0098,0.0121,0.0102,0.0113};
  double out[N] = {5.2538,5.2433,5.2599,5.2512,5.2398,5.2934,5.2799,5.2919};
  double outE[N] = {0.0121,0.0116,0.0306,0.0122,0.0099,0.0119,0.0095,0.0127};
  double x[N] = {1,2,3,4,5,6,7,8};
  double xE[N] = {0,0,0,0,0,0,0,0};
  TGraphErrors *gr = new TGraphErrors();
  for(int i=0;i<N;++i){
    gr->SetPoint(i, x[i], (in[i]-out[i])/5.25);
    gr->SetPointError(i,xE[i],sqrt(inE[i]*inE[i]+outE[i]*outE[i])/5.25);
  }
  gr->SetMarkerStyle(8);
  gr->SetTitle("In-Out Difference for CREX Moller Measurements");
  gr->GetXaxis()->SetTitle("Measurement Number");
  gr->GetYaxis()->SetTitle("#frac{|A_{in}|-|A_{out}|}{<A>}");
  gr->Draw("ap");
  gr->Fit("pol0");
  
}
