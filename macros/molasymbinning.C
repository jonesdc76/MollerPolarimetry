{
  long  N = 1e5*512;
  TRandom3 *r = new TRandom3(0);
  double RATE = 50, asym = 0.053, rate = RATE;
  TGraph *gr = new TGraph();
  TH1D *h1 = new TH1D("h1","h1",200,0.995,1.005);
  int s = 0;
  // while(rate < RATE){
  //   rate += 100;
  //   cout<<"Rate "<<rate<<endl;
  //   h = new TH1D("h","h",500,-1.15, 1.2);
  //   long n=0;
  //   for(int i=0;i<N;++i){
  //     double np = r->Gaus(rate*(1+asym),sqrt(rate*(1+asym)));
  //     double nn = r->Gaus(rate*(1-asym),sqrt(rate*(1-asym)));
  //     np = (np - int(np) > 0.5 ? double((int)np)+1 : double((int)np));
  //     nn = (nn - int(nn) > 0.5 ? double((int)nn)+1 : double((int)nn));
  //     n += nn+np;
  //     double pol = (double)(np-nn)/(double)(np+nn);
  //     h->Fill(pol);
  //     if( n > N ) break;
  //   } 
  //   gr->SetPoint(s, rate*2., h->GetMean()/asym);
  //   h1->Fill(h->GetMean()/asym);
  //   cout<<"Overflows: "<<h->GetBinContent(0)+h->GetBinContent(h->GetNbinsX()+1)<<endl;
  //   ++s;
  // }
  const int len = 10;
  double numer[len], denom[len];
  int npa[len] = {1,2,4,8,16,32,64,128,256,512};
  TH1D *h[len];
  double NperA = 20;

  for(int s=0;s<len;++s){
    numer[s] = 0;
    denom[s] = 0;
    h[s] = new TH1D(Form("h%i",s),Form("h%i",s),500,-1.15, 1.2);
  }
  long n=0;
  long i = 0;
  while(i < N){
    if(i%10000000==0)cout<<i<<endl;
    double np = r->Gaus(rate*(1+asym),sqrt(rate*(1+asym)));
    double nn = r->Gaus(rate*(1-asym),sqrt(rate*(1-asym)));
    np = (np - int(np) > 0.5 ? double((int)np)+1 : double((int)np));
    nn = (nn - int(nn) > 0.5 ? double((int)nn)+1 : double((int)nn));
    n += nn+np;
    for(int s=0;s<len;++s){
      numer[s] += np - nn;
      denom[s] += np + nn;
      if(i%npa[s] == 0){
	//if(s==0)cout<<numer[s]<<" "<<denom[s]<<endl;
	h[s]->Fill(numer[s]/denom[s]);
	numer[s] = 0; denom[s] = 0;
      }
    }
    ++i;
  } 
  for(int s=0;s<len;++s){
    gr->SetPoint(s, rate*npa[s], h[s]->GetMean()/asym);
    cout<<"Overflows: "<<h[s]->GetBinContent(0)+h[s]->GetBinContent(h[s]->GetNbinsX()+1)<<endl;
  }
  TCanvas *c = new TCanvas("c","c",0,0,1200,700);c->SetLogx();
  //  c->Divide(2,1);
  //  c->cd(1);
  //  c->cd(2);
  gr->SetMarkerStyle(8);
  gr->Draw("ap");
  gr->SetTitle("Simulated Ratio of Measured Asym/Input Asym vs N");
  gr->GetXaxis()->SetTitle("Number of Events per Asymmetry");
  gr->GetYaxis()->SetTitle("Measured Asym/Input Asym");

  


 
}
