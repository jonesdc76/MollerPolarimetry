{
  // const int N = 7;
  // double y4[N]={89.37, 89.52, -87.82, 89.12, (-88.10-87.78)/2.0, 89.30, -88.10};
  // double y4e[N]={0.29,  0.28,   0.28,  0.27,       0.27/sqrt(2),  0.29,   0.27};

  // double y10[N]={90.26, 89.92, -88.89, (89.95/0.32/0.32+90.77/0.28/0.28)/(1/0.32/0.32+1/0.28/0.28), -88.95, 90.50, -88.99};
  // double y10e[N]={0.23, 0.24,    0.27, 1/sqrt(1/(0.32*0.32)+1/(0.28*0.28)),                           0.24,  0.28,   0.27};
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);
  const int N = 9;
  double y4[N]={ 89.37, 89.52, -87.82,  89.12,  88.89, -87.78,  89.64, 89.30, -88.10};
  double y4e[N]={ 0.29,  0.28,   0.28,   0.27,   0.34,   0.27,   0.44,  0.29,   0.27};

  double y10[N]={90.26, 89.92, -88.89,  90.77,  89.95, -88.95, 90.08, 90.50, -88.99};
  double y10e[N]={0.23,  0.24,    0.27,  0.28,   0.32,   0.24,  0.23,  0.28,  0.27};
  //                                    74/71   73/66   63/70  62/98
  TGraphErrors *gr = new TGraphErrors();
  for(int i=0;i<N;++i){
    double x=i*1.0+1.0, y=y10[N-i-1]/y4[N-i-1];
    gr->SetPoint(i,x,y);
    double ye = y*sqrt(pow(y10e[N-i-1]/y10[N-i-1],2)+pow(y4e[N-i-1]/y4[N-i-1],2));
    gr->SetPointError(i,0,ye);
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
  TPaveText *ptt = new TPaveText(0.64,0.74,0.945,0.94,"ndc");
  ptt->SetBorderSize(0);
  ptt->SetFillColor(0);
  ptt->AddText(Form("#chi^{2}/NDF                    %0.2f/%i",f->GetChisquare(),f->GetNDF()));
  ptt->AddText(Form("P-value                   %0.4f",f->GetProb()));
  ptt->AddText(Form("Constant           %0.4f#pm%0.4f",f->GetParameter(0), f->GetParError(0)));
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  gr->Draw("ap");
  ptt->SetTextAlign(11);
  ptt->Draw("l");
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetTitleSize(0.04);
  //gr->GetYaxis()->SetTitleOffset(1.1);
  gr->GetXaxis()->SetTitle("Measurement Number");
  gr->GetYaxis()->SetTitle("Asymmetry Ratio (10 #mum)/(4 #mum)");
  c->SaveAs("ratio10to4foil.pdf");
}
