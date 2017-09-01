{
  //This summarizes the g-factor values for Fe in 
  //(1)Bagguley "Ferromagnetic resonance in colloidal suspension" 1953
  //(2)Rodbell "Ferromagnetic resonance of iron whisker crystals" 1959
  //(2)Table 2 Meyer and Asch "Experimental g-prime and g Values of Fe, Co, Ni 
  //and Their Alloys" 1961
  //(3) Frait "The g-factor in pure polycrystalline iron" 1977
  //(4) Frait "THE g-FACTOR AND SURFACE MAGNETIZATION OF PURE IRON ALONG [100] AND [111] DIRECTIONS" 1971
  const int N=8;
  int color[N] = {kGray+2,kBlack,kGreen+3,kRed,kRed,kRed,kBlue,kRed};
  int style[N] = {34,21,8,4,4,4,33,34};
  TCanvas *c = new TCanvas("c","c",0,0,700,500);
  TString auth[N] = {"Bagguley 1953","Barlow & Standley 1956", "Rodbell 1959","Meyer & Asch 1961","Meyer & Asch 1961","Meyer & Asch 1961", "Frait & Gemperle 1971","Frait 1971"};
  double x[N+2] = {2.16, 2.10,2.05,2.092,2.090,2.093,2.089,2.088}, xe[N+2]={0.0216,0.020,0.010, 0.015,0.015,0.015,0.008,0.008};
  double y[N+2]={1,2,3,4,5,6,7,8},ye[N+2]={0,0,0,0,0,0,0,0};
  TLegend *leg = new TLegend(0.6, 0.89, 0.9, 0.5);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetShadowColor(0);
  TGraphErrors *gr[N];
  TMultiGraph *mg = new TMultiGraph();
  for(int n=0;n<N;++n){
    gr[n] = new TGraphErrors(1,&(x[n]),&(y[n]),&(xe[n]),&(ye[n]));
    gr[n]->SetMarkerColor(color[n]);
    gr[n]->SetMarkerStyle(style[n]);
    mg->Add(gr[n]);
    if(n!=4&&n!=5)
      leg->AddEntry(gr[n], auth[n].Data(),"p");
  }

  mg->Draw("ap");
  mg->SetTitle("Compiled Measurements of Spectroscopic g-factor for Fe");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("g\'");
  mg->GetXaxis()->SetLimits(2.03,2.19);
  //  mg->GetXaxis()->SetRangeUser(2.03,2.19);
  leg->Draw();
  double fitval = 2.08579;
  TLine tl = TLine(fitval,mg->GetYaxis()->GetXmin(),fitval,mg->GetYaxis()->GetXmax());
  tl.SetLineColor(kBlack);
  tl.SetLineWidth(2);
  tl.Draw();
  c->SaveAs("gfactor_world_data.pdf");
}
FARS_RPC_GAS_101	pressure_sf6	High!!!	30.96	2017-09-01 03:39:31
FARS_RPC_GAS_102	HV_interlock	Low!!!	0.00	2017-09-01 03:39:31
FARS_RPC_GAS_102	flowrate_status	Low!!!	0.00	2017-09-01 03:39:31
FARS_RPC_GAS_101	pressure_argon	High!!!	30.00	2017-09-01 03:02:11
FARS_RPC_GAS_101	pressure_r134a	High!!!	30.14	2017-09-01 02:16:01
FARS_RPC_GAS_101	flow_rate_r134a	Low!!!	523.60	2017-09-01 01:58:21
FARS_RPC_GAS_101	flow_rate_isobutane	Low!!!	69.76	2017-09-01 01:57:50
FARS_RPC_GAS_101	flow_rate_argon	Low!!!	1140.00	2017-09-01 01:57:41
FARS_RPC_GAS_101	flow_rate_sf6	Low!!!	8.87	2017-09-01 01:57:41
FARS_RPC_GAS_101	pressure_isobutane	Low!!!	19.92	2017-09-01 00:41:01
FARS_RPC_GAS_102	HAD_2_OUT	High!!!	63.70	2017-08-31 15:49:40
