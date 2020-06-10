{
  TH1D *h = new TH1D("h","h",100,-1,1);
  TCanvas *c = new TCanvas("c","c",0,0,1700,700);
  c->Divide(2,1);
  bool linear = 0;
  c->cd(1);
  TLegend *leg = new TLegend(0.8,0.1,0.99,0.9);
   const UInt_t Number = 3;
   const bool Eugene = 0;
   const double R = 5.08 - Eugene*0.08, R0 = 3.7;
   Double_t Red[Number]    = { 1.00, 0.00, 0.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   const int n = 22;
   //double param[10] = {0,0,-0.00035,-0.00036,0.00116,-0.01441,0.00012,0.00003,0.00007,0.00056};
   double param[10] = {0,0,0,0,0,-0.01441,0,0,0,0};
   for(int i=2; i<10;++i)param[i] *= pow(R/R0,i-1);
   TMultiGraph *mg = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();
   TGraph *gr[2*n];
   double cur[n];
   TString s = "%lg";
   TF1 *f[n];
   for(int i = 0; i<n; ++i){
     cur[i] = 300 - i*30;
     gr[i] = new TGraph("QO1H01.dat", Form("%s%s",s.Data(),"%lg"));
     mg->Add(gr[i]);
     int col = TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,n) + i;
     gr[i]->SetMarkerColor(col);
     gr[i]->SetMarkerStyle(21+(i%15));
     gr[i]->SetMarkerSize(0.9);
     s += "%*lg";
     if(linear){
       f[i] = new TF1(Form("f%i",i),"pol1",0,0);
     }else{
       // f[i] = new TF1(Form("f%i",i),Form(
       // 		      "[0]*(x-[1])*(1+[2]*pow((x-[1])/%f,1)+[3]*pow((x-[1])/%f,2)"
       // 		      "+[4]*pow((x-[1])/%f,3)+[5]*pow((x-[1])/%f,4)"
       // 		      "+[6]*pow((x-[1])/%f,5)+[7]*pow((x-[1])/%f,6)"
       // 		      "+[8]*pow((x-[1])/%f,7)+[9]*pow((x-[1])/%f,8))",R,R,R,R,R,R,R,R),-5,5);
       f[i] = new TF1(Form("f%i",i),Form("[0]*(x-[1])*(1+[2]*pow((x-[1])/%f,4))",R),-5,5);
       f[i]->SetParameter(2,-0.048);
       if(Eugene)
	 f[i]->FixParameter(2,-0.048);
       // f[i]->SetParameters(param);
       // f[i]->SetParameter(5,-0.12);
       // f[i]->FixParameter(2,param[2]);
       // f[i]->FixParameter(3,param[3]);
       // f[i]->FixParameter(4,param[4]);
       // // f[i]->FixParameter(5,param[5]);
       // f[i]->FixParameter(6,param[6]);
       // f[i]->FixParameter(7,param[7]);
       // f[i]->FixParameter(8,param[8]);
       // f[i]->FixParameter(9,param[9]);
     }
     f[i]->SetLineColor(col);
     f[i]->SetLineStyle(1+(i%10));
     gr[i]->SetLineColor(col);
     gr[i]->SetLineStyle(1+(i%10));
     double x, y;
     gr[i]->GetPoint(i,x,y);
     TFitResultPtr ftr = gr[i]->Fit(f[i],"S");
     cout<<cur[i]<<" "<<int(ftr)<<endl;
     while(int(gr[i]->Fit(f[i],"S"))!=0){
       cout<<"Trying again "<<endl;
     }

     if(linear)
        h->Fill(-f[i]->GetParameter(0)/f[i]->GetParameter(1));
     else{
       h->Fill(f[i]->GetParameter(1));       
       //h->Fill(f[i]->GetParameter(5));
     }
     gr[i+n] = new TGraph();
     for(int j=0;j<gr[i]->GetN();++j){
       double x, y;
       gr[i]->GetPoint(j, x,y);
       gr[i+n]->SetPoint(j, x, 100*(y-f[i]->Eval(x))/y);
     }
     gr[i+n]->RemovePoint(3);
     leg->AddEntry(gr[i],Form("%i Amps",(int)cur[i]),"lp");
     mg2->Add(gr[i+n]);
     gr[i+n]->SetLineColor(col);
     gr[i+n]->SetLineStyle(1+(i%10));
     gr[i+n]->SetMarkerColor(col);
     gr[i+n]->SetMarkerStyle(21+(i%15));
     gr[i+n]->SetMarkerSize(0.9);
   }
   
   mg->Draw("ap");
   mg->SetTitle("QO1H01 (Q1) Bdl vs X-position");
   mg->GetXaxis()->SetTitle("X-position (cm)");
   mg->GetYaxis()->SetTitle("Bdl (Gauss cm)");
   mg->GetYaxis()->SetTitleOffset(1.35);
   mg->GetXaxis()->SetLimits(-6,6);
   mg->GetXaxis()->SetRangeUser(-4.5, 5.5);
   leg->Draw();
   c->cd(2);
   mg2->Draw("alp");
   mg2->SetTitle("QO1H01 (Q1) Residual of Linear Fit");
   mg2->GetXaxis()->SetTitle("X-position (cm)");
   mg2->GetYaxis()->SetTitle("Residual (%)");
   mg2->GetYaxis()->SetTitleOffset(0.92);
   mg2->GetXaxis()->SetLimits(-6,6);
   mg2->GetXaxis()->SetRangeUser(-4.5, 5.5);
   leg->Draw();
   gPad->SetGrid();
   new TCanvas;
   h->Draw();

   //Build latex table
   printf("Current &Offset & Dodecapole \\\\ \n");
   printf("(A)     &(mm)   &Coeff (frac)\\\\ \\hline \n");
   double p1 = 0, p2 = 0;
   for(int i=0;i<n;++i){
     printf("%7.1f & %5.2f & %10.3f \\\\ \\hline\n",cur[i],f[i]->GetParameter(1)*10,f[i]->GetParameter(2));
     p1 += f[i]->GetParameter(1)*10/double(n);
     p2 += f[i]->GetParameter(2)/double(n);
   }
   printf("\\hline\n");
   printf("Average & %5.2f & %10.3f \\\\ \\hline \n",p1,p2);
   if(Eugene) c->SaveAs("EugeneCompNewQ1.pdf");
   else c->SaveAs("QO1H01_BdlvsX.pdf");
}
