{
  //This is a combination of the tabulated values in Table 1 of these two refs
  //(1)G. G. Scott "Gryomagnetic Ratio Experiments" 1962 and 
  //(2)Meyer and Asch "Experimental g-prime and g Values of Fe, 
  //Co, Ni and Their Alloys" 1961.
  //The error on Scott's gprime value 1.835+/-0.002 was inflated to +/-0.008
  //due to the the effect of impurities. His Ni sample had ~0.2% impurities and 
  //given the potentially large effect of these (sample of Barnett and Kenney
  //with 1.4% impurities had gprime that was 3% larger than values measured by 
  //Scott and Meyer) I assigned an additional systematic error of 0.4% in quad-
  //rature with Scott's 0.002. The error in Meyer's 1957 measurement was NOT 
  //inflated from +/-0.010 even though there were 0.1% impurities due to its
  //magnetic properties (Curie T and saturation magnetization) not being 
  //apparently affected by the impurities see Table 2 in (1)

  const int N=4;
  int color[N] = {kBlack,kRed,kBlue,kBlack};
  int style[N] = {34,21,8,34};
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TString auth[N] = { "Proposed Error", "Scott 1952-1960", "Meyer 1957", "Meyer 1958", };
  double x[N] = {1.835, 1.835, 1.852, 1.845}, xe[N]={0.006, 0.002, 0.01, 0.008};
  double y[N],ye[N];
  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.15);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  TGraphErrors *gr[N];
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<N;++i){
    y[i] = (i==0? 1: i);
    ye[i] = 0;
    gr[i] = new TGraphErrors(1,&(x[i]),&(y[i]),&(xe[i]),&(ye[i]));
    gr[i]->SetMarkerColor(color[i]);
    gr[i]->SetLineWidth(2);
    if(i==0){
      gr[i]->SetLineColor(color[i]);
      gr[i]->SetLineStyle(2);
    }
    gr[i]->SetMarkerStyle(style[i]);
    mg->Add(gr[i]);
    if(i==0)
      leg->AddEntry(gr[i], auth[i].Data(),"l");
    else
      leg->AddEntry(gr[i], auth[i].Data(),"p");
   }
  mg->Draw("ap");
  mg->SetTitle("Compiled Measurements of g\' for Ni");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("g\'");
  //  mg->GetXaxis()->SetLimits(1.82,1.88);
  leg->Draw();
  double fitval = 1.841;//1.84111+/-0.0043 chsq 2.459/2 with proposed error
  //double fitval = 1.836;//1.83618//+/-0.0019 chsq 4.066/2 with Scott's error
  TLine tl = TLine(fitval,mg->GetYaxis()->GetXmin(),fitval,mg->GetYaxis()->GetXmax());
  tl.SetLineColor(kBlack);
  tl.SetLineWidth(2);
  tl.Draw();
  c->SaveAs("gprime_world_data_Ni.pdf");
}
