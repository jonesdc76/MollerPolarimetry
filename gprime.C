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
  const int N=8;
  int color[N] = {kGray+2,kBlack,kGreen+3,kRed,kBlue,kRed,kViolet,kOrange+7};
  int style[N] = {34,21,8,4,33,34,26,21};
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TString auth[N] = {"Barnett 1944","Meyer 1951", "Scott 1951", "Barnett & Kenny 1952", "Scott 1955","Meyer & Brown 1957", "Scott (cylinder) 1960","Scott (ellipsoid) 1960"};
  double x[N] = {1.937,1.936,1.927,1.929,1.919,1.929,1.917,1.919}, xe[N]={0.006,0.008,0.004,0.006,0.006,0.008,0.002,0.002};
  double y[N],ye[N];
  TLegend *leg = new TLegend(0.6, 0.89, 0.9, 0.5);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  TGraphErrors *gr[N];
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<N;++i){
    y[i] = i+1;
    ye[i] = 0;
    gr[i] = new TGraphErrors(1,&(x[i]),&(y[i]),&(xe[i]),&(ye[i]));
    gr[i]->SetMarkerColor(color[i]);
    gr[i]->SetMarkerStyle(style[i]);
    mg->Add(gr[i]);
    leg->AddEntry(gr[i], auth[i].Data(),"p");
  }
  mg->Draw("ap");
  mg->SetTitle("Compiled Measurements of g\' for Fe");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("g\'");
  mg->GetXaxis()->SetLimits(1.91,1.96);
  mg->GetXaxis()->SetRangeUser(1.91,1.96);
  leg->Draw();
  double fitval = 1.9208;
  TLine tl = TLine(fitval,mg->GetYaxis()->GetXmin(),fitval,mg->GetYaxis()->GetXmax());
  tl.SetLineColor(kBlack);
  tl.SetLineWidth(2);
  tl.Draw();
  c->SaveAs("gprime_world_data_Fe.pdf");
}
