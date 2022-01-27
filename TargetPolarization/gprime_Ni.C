{
  //This is a combination of the tabulated values in Table 1 of these two refs
  //(1)G. G. Scott "Gryomagnetic Ratio Experiments" 1962 and 
  //(2)Meyer and Asch "Experimental g-prime and g Values of Fe, 
  //Co, Ni and Their Alloys" 1961.

  const int N=3;
  int color[N+1] = {kRed,kBlue,kBlack,kBlack};
  int style[N+1] = {21,8,34,1};
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TString auth[N] = { "Scott 1952-1960", "Meyer 1957", "Meyer 1958"};
  double x[N+1] = {1.835, 1.852, 1.845,1.861}, xe[N+1]={0.002, 0.01, 0.008,0.002*1.836};
  double y[N+1],ye[N+1];
  TPaveText *pt = new TPaveText(0.68,0.63,0.765,0.68,"ndc");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  pt->AddText("0.4%");
  TLegend *leg = new TLegend(0.6, 0.89, 0.89, 0.7);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  TGraphErrors *gr[N+1];
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<=N;++i){
    y[i] = i==N? 2.57:N-i;
    ye[i] = 0;
    gr[i] = new TGraphErrors(1,&(x[i]),&(y[i]),&(xe[i]),&(ye[i]));
    gr[i]->SetMarkerColor(color[i]);
    //gr[i]->SetLineColor(color[i]);
    gr[i]->SetMarkerStyle(style[i]);
    mg->Add(gr[i]);
    if(i<N)
      leg->AddEntry(gr[i], auth[i].Data(),"pp");
   else
      gr[i]->SetLineWidth(2);
  }
  mg->GetYaxis()->SetLimits(0.5,3.5);
  mg->GetYaxis()->SetRangeUser(0.5,3.5);
  mg->SetTitle("Compiled Measurements of g\' for Ni");
  mg->Draw("ap");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("g\'");
  mg->GetXaxis()->SetLimits(1.83,1.87);
  leg->Draw();
  pt->Draw();
  double fitval = 1.8362;//1.83618//+/-0.0019 chsq 4.066/2 with Scott's error
  TLine tl = TLine(fitval,mg->GetYaxis()->GetXmin(),fitval,mg->GetYaxis()->GetXmax());
  tl.SetLineColor(kBlack);
  tl.SetLineWidth(2);
  tl.Draw();
  c->SaveAs("gprime_world_data_Ni.pdf");
  c->SaveAs("../nim/figures/gprime_world_data_Ni.pdf");
  for(int i=0;i<N;++i)xe[i]*=1.72;
  //xe[1]*=3;
  TCanvas *c1 = new TCanvas("c1","c1",700,500,600,400);
  TGraphErrors *grx = new TGraphErrors(N,y,x,ye,xe);
  gStyle->SetOptFit(1111);
  grx->SetMarkerStyle(8);
  grx->Draw("ap");
  grx->Fit("pol0");

}
