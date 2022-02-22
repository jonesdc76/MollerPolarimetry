{
  // const int N = 7;
  // double y4[N]={89.37, 89.52, -87.82, 89.12, (-88.10-87.78)/2.0, 89.30, -88.10};
  // double y4e[N]={0.29,  0.28,   0.28,  0.27,       0.27/sqrt(2),  0.29,   0.27};

  // double y10[N]={90.26, 89.92, -88.89, (89.95/0.32/0.32+90.77/0.28/0.28)/(1/0.32/0.32+1/0.28/0.28), -88.95, 90.50, -88.99};
  // double y10e[N]={0.23, 0.24,    0.27, 1/sqrt(1/(0.32*0.32)+1/(0.28*0.28)),                           0.24,  0.28,   0.27};
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);
  const int N = 8;
  double a65 = -5.4515, a65e = 0.0170, a70 = -5.4319, a70e = 0.017, a63 = -5.5046, a63e = 0.015,
    a62 = 5.5746, a62e = 0.0144, a66 = 5.5007, a66e = 0.0213, a73 = 5.5666, a73e = 0.0199, a74 = 5.6169, a74e = 0.0176, a71 = 5.5147, a71e = 0.0167, a98 = 5.5469, a98e = 0.0275;
  double a1d = (a65/pow(a65e,2)+a70/pow(a70e,2))/(pow(a65e,-2)+pow(a70e,-2));
  double a1de = sqrt(1/(pow(a65e,-2)+pow(a70e,-2)));
  double a2n = ( a73/pow(a73e,2)+a74/pow(a74e,2) )/(pow(a73e,-2)+pow(a74e,-2)) ,
    a2ne = 1/sqrt(pow(a73e,-2)+pow(a74e,-2));
  cout<<a2n<<" "<<a2ne<<endl;

  double a2d = ( a71/pow(a71e,2)+a98/pow(a98e,2) ) / (pow(a71e,-2)+pow(a98e,-2)) ,
        a2de = 1/sqrt(pow(a71e,-2)+pow(a98e,-2));
  cout<<a2d<<" "<<a2de<<endl;
  
  double y4[N] = {  5.5304,  -5.4345, 5.5394,  a1d,   a2d,    a66, -5.4517,  5.5260};
  double y4e[N] = { 0.0176,   0.0170, 0.0174, a1de,  a2de,   a66e,  0.0170,  0.0182};
			   	     	     	  		            
  double y10[N] = { 5.5853,  -5.5006, 5.5643,  a63,   a2n,    a62, -5.5067,  5.6000};
  double y10e[N] = {0.0145,   0.0166, 0.0151, a63e,  a2ne,   a62e,  0.0164,  0.0172};
  
  TGraphErrors *gr = new TGraphErrors();
  for(int i=0;i<N;++i){
    double x=i*1.0+1.0, y=y10[N-i-1]/y4[N-i-1];
    gr->SetPoint(i, x, y);
    double ye = y*sqrt(pow(y10e[N-i-1]/y10[N-i-1],2)+pow(y4e[N-i-1]/y4[N-i-1],2));
    gr->SetPointError(i, 0, ye);
  }
  gr->SetMarkerStyle(8);
  gStyle->SetOptFit(0);
  gStyle->SetFitFormat("6.4g");
  TF1 *f = new TF1("f","pol0",0,10);
  f->SetLineWidth(3);
  f->SetLineColor(kBlack);
  f->SetLineStyle(10);
  f->SetParNames("Constant");
  gr->Fit(f);
  TPaveText *ptt = new TPaveText(0.62,0.74,0.945,0.94,"ndc");
  ptt->SetBorderSize(0);
  ptt->SetFillColor(0);
  ptt->AddText(Form("#chi^{2}/NDF                    %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  ptt->AddText(Form("P-value                   %0.4f",f->GetProb()));
  ptt->AddText(Form("Constant           %0.4f #pm %0.4f",f->GetParameter(0), f->GetParError(0)));
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  gr->Draw("ap");
  ptt->SetTextAlign(11);
  ptt->Draw("l");
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetRangeUser(0.999,1.023);
  //gr->GetYaxis()->SetTitleOffset(1.1);
  gr->GetXaxis()->SetTitle("Measurement Number");
  gr->GetYaxis()->SetTitle("Asymmetry Ratio (10 #mum)/(4 #mum)");
  gPad->Update();
  c->SaveAs("ratio10to4foil.pdf");
}
