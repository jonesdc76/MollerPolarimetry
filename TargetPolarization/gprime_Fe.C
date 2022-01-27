{
  //This is a combination of the tabulated values in Table 1 of these two refs
  //(1)G. G. Scott "Gryomagnetic Ratio Experiments" 1962 and 
  //(2)Meyer and Asch "Experimental g-prime and g Values of Fe, 
  //Co, Ni and Their Alloys" 1961
  //There were two inconsistencies between the two Tables that were resolved
  //as follows: I. Table 1 of (2) has Barnett 1941 rho e/mc=1.035 whereas 
  //Barnett actually had the following 3 values for different Fe samples: 
  //1.032, 1.032 and 1.034. (2) used 1.032 for Barnett 1944. I use the straight
  //average of 1.0327 or gprime= 1.937+/-0.006 where error is taken from (1). 
  //II. (2) has gprime = 1.932 for Meyer and Brown 1957 which is inconsistent 
  //with what Meyer himself gives in (2). For this I use the value from (2) 
  //along with the error on this value from (1) gprime=1.929+/-0.008
  const int N=9;
  int color[N] = {kGray+2,kBlack,kGreen+3,kRed,kBlue,kRed,kViolet,kOrange+7,kBlack};
  int style[N] = {34,21,8,4,33,34,26,21};
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TString auth[N] = {"Barnett 1944","Meyer 1951", "Scott 1951", "Barnett & Kenny 1952", "Scott 1955","Meyer & Brown 1957", "Scott (cylinder) 1960","Scott (ellipsoid) 1960"};
  double x[N] = {1.937,1.936,1.927,1.929,1.919,1.929,1.917,1.919,1.95}, xe[N]={0.006,0.008,0.004,0.006,0.006,0.008,0.002,0.002,0.002*1.9208};
  double y[N],ye[N];
  TLegend *leg = new TLegend(0.6, 0.3, 0.89, 0.7);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  TPaveText *pt = new TPaveText(0.70,0.22,0.78,0.265,"ndc");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  pt->AddText("0.4%");
  TGraphErrors *gr[N];
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<N;++i){
    y[i] = i==N-1 ? 2:N-i-1;
    ye[i] = 0;
    gr[i] = new TGraphErrors(1,&(x[i]),&(y[i]),&(xe[i]),&(ye[i]));
    gr[i]->SetMarkerColor(color[i]);
    gr[i]->SetMarkerStyle(style[i]);
    mg->Add(gr[i]);
    if(i!=N-1)
      leg->AddEntry(gr[i], auth[i].Data(),"p");
    else
      gr[i]->SetLineWidth(2);

  }
  
  mg->Draw("ap");
  mg->SetTitle("Compiled Measurements of g\' for Fe");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("g\'");
  mg->GetXaxis()->SetLimits(1.91,1.96);
  mg->GetXaxis()->SetRangeUser(1.91,1.96);
  leg->Draw();
  pt->Draw();
  gPad->Update();
  double fitval = 1.9208;
  TLine tl = TLine(fitval,mg->GetYaxis()->GetXmin(),fitval,mg->GetYaxis()->GetXmax());
  tl.SetLineColor(kBlack);
  tl.SetLineWidth(2);
  tl.Draw();
  // TGraphErrors *gt = new TGraphErrors(8,y,x,ye,xe);
  // gt->Draw("ap");
  // gStyle->SetOptFit(1111);
  // gt->Fit("pol0");
  c->SaveAs("gprime_world_data_Fe.pdf");
  c->SaveAs("../nim/figures/gprime_world_data_Fe.pdf");
  for(int i=0;i<N;++i)xe[i]*=1.81;
  TCanvas *c1 = new TCanvas("c1","c1",700,500,600,400);
  TGraphErrors *grx = new TGraphErrors(N-1,y,x,ye,xe);
  gStyle->SetOptFit(1111);
  grx->SetMarkerStyle(8);
  grx->Draw("ap");
  TF1 *f = new TF1("f","pol0",0.9,9);
  grx->Fit("f","r");
  cout<<"Probability: "<<f->GetProb()<<endl;
  
}
