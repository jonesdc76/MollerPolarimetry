{
  const int nx = 8;
   const double* par1[nx] = {slope1,slope2,slope3,slope4,slope5,slope6,slope7,slope8};
   double err[nx] = {sloperr1,sloperr2,sloperr3,sloperr4,sloperr5,sloperr6,sloperr7,sloperr8};
   double num[nx] = {1,2,3,4,5,6,7,8};
   string date[nx] = {"April","End of May","End of June","Mid July","Beg of Aug","End of Sept","Beg of Oct","Beg of Oct"};
   Long64_t n = 8;
   double mean = TMath::Mean(n,par1);
   //  TH2F *h = new TH2F("h","variation of center position",nx,0.0,nx,nx,-0.0005,0.00045);
   //   TH1D *h = new TH1D("h","variation of center position",nx,1,nx);
   //   h->SetBit(TH1::kCanRebin);
   //   h->SetStats(0);
   //  double x[nx],y[nx];

   TGraphErrors* gr1 = new TGraphErrors(nx,num,par1,0,err);
   TF1 *mfunc = new TF1("mfunc","mean + [0]*x",0,8);

     for(int i=0;i<=nx;i++) {
         x[i]=i; 
        y[i]=par1[i];
       h->Fill(i,par1[i]);
         h->Fill(par1[i]);
         h->SetBinError(i,err[i]);
     }
  }

   gr1->SetMarkerColor(4);
   gr1->SetMarkerStyle(21);
   gr1->Draw("ALP");
   mfunc->Draw("SAME");
}
