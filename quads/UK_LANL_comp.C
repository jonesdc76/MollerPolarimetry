{
  //Patsy comparison between UK and LANL
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.46);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);
  TCanvas *cp = new TCanvas("cp","cp",0,0,700,500);
  TMultiGraph *mgP = new TMultiGraph();
  TGraph *grP = new TGraph("LANL_GLvsCur.dat","%lg%lg");
  TGraph *grPuk = new TGraph("UK_GLvsCur.dat","%lg%lg");
  gStyle->SetOptFit(1111);
  for(int i=0;i<grP->GetN();++i){
    double x, y;
    grP->GetPoint(i,x,y);
    grP->SetPoint(i,x/300.,y);
  }
  for(int i=0;i<grPuk->GetN();++i){
    double x, y;
    grPuk->GetPoint(i,x,y);
    grPuk->SetPoint(i,x/300.,y/10000.);
  }
  TF1 *fp = new TF1("fp","pol5",-300,300);
  fp->SetLineColor(kBlue);
  grPuk->Fit(fp);
  grP->SetMarkerStyle(21);
  grPuk->SetMarkerColor(kBlue);
  grPuk->SetMarkerStyle(8);
  mgP->Add(grPuk);
  mgP->Add(grP);
  TGraph *grPdp = new TGraph();
  TGraph *grPd = new TGraph();
  for(int i=0;i<grP->GetN()-1;++i){
    double x, y;
    grP->GetPoint(i,x,y);
    grPdp->SetPoint(i,x,100*(y-fp->Eval(x))/y);
    grPd->SetPoint(i,x,100*(y-fp->Eval(x)));
    cout<<x<<" "<<y<<endl;
  }
  grPd->SetMarkerStyle(8);
  grPd->SetMarkerColor(kRed);
  grPdp->SetMarkerStyle(23);
  grPdp->SetMarkerColor(kGreen+2);
  mgP->Add(grPd);
  mgP->Add(grPdp);
  mgP->Draw("ap");
  mgP->SetTitle("Patsy Gradient-Length versus Current");
  mgP->GetXaxis()->SetTitle("Current/300 A");
  mgP->GetYaxis()->SetTitle("GL (T)");
  TLegend *tl = new TLegend(0.15,0.7,0.32,0.87);
  tl->AddEntry(grP,"LANL","p");
  tl->AddEntry(grPuk,"UK","lp");
  tl->AddEntry(grPdp, "Diff (#times100 T)","p");
  tl->AddEntry(grPd, "Diff(%)","p");
  gPad->SetGrid();
  tl->Draw();
  gPad->Update();
  double l = mgP->GetYaxis()->GetXmin();
  double h = mgP->GetYaxis()->GetXmax();
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),l,h,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetTitleColor(kRed);
  axis->SetTitle("LANL-UKPol5 (%)");
  axis->Draw();
  cp->SaveAs("Pat_LANL_UK_GLcomp.pdf");

  //Tessa comparison between UK and LANL
  TCanvas *ct = new TCanvas("ct","ct",0,0,700,500);
  TMultiGraph *mgT = new TMultiGraph();
  TGraph *grT = new TGraph("LANL_GLvsCur.dat","%*lg%*lg%lg%lg");
  TGraph *grTuk = new TGraph("UK_GLvsCur.dat","%lg%*lg%lg");
  for(int i=0;i<grT->GetN();++i){
    double x, y;
    grT->GetPoint(i,x,y);
    grT->SetPoint(i,x/300.,y);
  }
  for(int i=0;i<grTuk->GetN();++i){
    double x, y;
    grTuk->GetPoint(i,x,y);
    grTuk->SetPoint(i,x/300.,y/10000.);
  }
  TF1 *ft = new TF1("ft","pol5",-300,300);
  ft->SetLineColor(kBlue);
  grTuk->Fit(ft);
  grT->SetMarkerStyle(21);
  grTuk->SetMarkerColor(kBlue);
  grTuk->SetMarkerStyle(8);
  mgT->Add(grTuk);
  mgT->Add(grT);
  TGraph *grTd = new TGraph();
  TGraph *grTdp = new TGraph();
  for(int i=0;i<grT->GetN()-1;++i){
    double x, y;
    grT->GetPoint(i,x,y);
    grTd->SetPoint(i,x,100*(y-ft->Eval(x)));
    grTdp->SetPoint(i,x,100*(y-ft->Eval(x))/y);
  }
  grTdp->SetMarkerStyle(8);
  grTdp->SetMarkerColor(kRed);
  mgT->Add(grTdp);
  grTd->SetMarkerStyle(22);
  grTd->SetMarkerColor(kGreen+2);
  mgT->Add(grTd);
  mgT->Draw("ap");
  mgT->SetTitle("Tessa Gradient-Length versus Current");
  mgT->GetXaxis()->SetTitle("Current/300 A");
  mgT->GetYaxis()->SetTitle("GL (T)");
  gPad->SetGrid();
  tl->Draw();
  gPad->Update();
  l = mgT->GetYaxis()->GetXmin();
  h = mgT->GetYaxis()->GetXmax();
  TGaxis *axisT = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),l,h,510,"+L");
  axisT->SetLineColor(kRed);
  axisT->SetLabelColor(kRed);
  axisT->SetTitleColor(kRed);
  axisT->SetTitle("LANL-UKPol5 (%)");
  axisT->Draw();
  ct->SaveAs("Tes_LANL_UK_GLcomp.pdf");

  //Felicia comparison between UK and LANL
  TCanvas *cf = new TCanvas("cf","cf",0,0,700,500);
  TMultiGraph *mgF = new TMultiGraph();
  TGraph *grF = new TGraph("LANL_GLvsCur.dat","%*lg%*lg%*lg%*lg%lg%lg");
  TGraph *grFuk = new TGraph("UK_GLvsCur.dat","%lg%*lg%*lg%lg");
  for(int i=0;i<grF->GetN();++i){
    double x, y;
    grF->GetPoint(i,x,y);
    grF->SetPoint(i,x/300.,y);
  }
  for(int i=0;i<grFuk->GetN();++i){
    double x, y;
    grFuk->GetPoint(i,x,y);
    grFuk->SetPoint(i,x/300.,y/10000.);
  }
  TF1 *ff = new TF1("ff","pol5",-300,300);
  ff->SetLineColor(kBlue);
  grFuk->Fit(ff);
  grF->SetMarkerStyle(21);
  grFuk->SetMarkerColor(kBlue);
  grFuk->SetMarkerStyle(8);
  mgF->Add(grFuk);
  mgF->Add(grF);
  TGraph *grFd = new TGraph();
  TGraph *grFdp = new TGraph();
  for(int i=0;i<grF->GetN()-1;++i){
    double x, y;
    grF->GetPoint(i,x,y);
    grFd->SetPoint(i,x,100*(y-ff->Eval(x)));
    grFdp->SetPoint(i,x,100*(y-ff->Eval(x))/y);
  }
  grFd->SetMarkerStyle(22);
  grFd->SetMarkerColor(kGreen+2);
  mgF->Add(grFd);
  grFdp->SetMarkerStyle(8);
  grFdp->SetMarkerColor(kRed);
  mgF->Add(grFdp);
  mgF->Draw("ap");
  mgF->SetTitle("Felicia Gradient-Length versus Current");
  mgF->GetXaxis()->SetTitle("Current/300 A");
  mgF->GetYaxis()->SetTitle("GL (T)");
  gPad->SetGrid();
  tl->Draw();
  gPad->Update();
  l = mgF->GetYaxis()->GetXmin();
  h = mgF->GetYaxis()->GetXmax();
  TGaxis *axisF = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),l,h,510,"+L");
  axisT->SetLineColor(kRed);
  axisT->SetLabelColor(kRed);
  axisT->SetTitleColor(kRed);
  axisT->SetTitle("LANL-UKPol5 (%)");
  axisT->Draw();
  cf->SaveAs("Fel_LANL_UK_GLcomp.pdf");
  
  TCanvas *cxt = new TCanvas("cxt","cxt",0,0,700,500);
  TMultiGraph *mgxt = new TMultiGraph();
  TGraph *grxtG4 = new TGraph();
  TGraph *grxtN = new TGraph();
  //  grTuk->Draw();
  TF1 *ftxG4 = new TF1("ftxG4",Form("0.000632446+5.15178*x-0.00262778*pow(x,2)"
				    "-0.107634*pow(x,3)+0.00209902*pow(x,4)"
				    "-0.640635*pow(x,5)"),-1,1);
  int n=0;
  for(int i=0;i<grTuk->GetN();++i){
    double x,y;
    grTuk->GetPoint(i,x,y);
    if(abs(x)<0.01)continue;
    
    grxtN->SetPoint(n,x,(y-ft->Eval(x))/y*100);
    grxtG4->SetPoint(n,x,(y-ftxG4->Eval(x))/y*100);
    ++n;
  }
  grxtG4->SetLineColor(kBlue);
  grxtN->SetLineColor(kGreen+2);
  grxtG4->SetMarkerColor(kBlue);
  grxtN->SetMarkerColor(kGreen+2);
  grxtG4->SetMarkerStyle(23);
  grxtN->SetMarkerStyle(8);
  mgxt->SetTitle("Tessa  Residual vs. Current (UK Data - Parametrization)");
  mgxt->GetXaxis()->SetTitle("Current/300 A");
  mgxt->GetYaxis()->SetTitle("Diff: UK Data - Pol5 (%)");
  gPad->SetGrid();
  mgxt->Add(grxtG4);
  mgxt->Add(grxtN);
  mgxt->Draw("acp");
  TLegend *tl2 = new TLegend(0.15,0.7,0.32,0.87);
  tl2->AddEntry(grxtG4, "G4 Pol5","lp");
  tl2->AddEntry(grxtN, "New Pol5","lp");
  tl2->Draw();
  cxt->SaveAs("Tessa_G4_UK_residual.pdf");
  
  TCanvas *cxf = new TCanvas("cxf","cxf",0,0,700,500);
  TMultiGraph *mgxf = new TMultiGraph();
  TGraph *grxfG4 = new TGraph();
  TGraph *grxfN = new TGraph();
  TF1 *ffxG4 = new TF1("ffxG4",Form("0.0001732+5.2119*x-0.000732518*pow(x,2)"
				    "-0.133423*pow(x,3)+0.000618402*pow(x,4)"
				    "-0.647082*pow(x,5)"),-1,1);
  n=0;
  for(int i=0;i<grFuk->GetN();++i){
    double x,y;
    grFuk->GetPoint(i,x,y);
    if(abs(x)<0.01)continue;
    
    grxfN->SetPoint(n,x,(y-ff->Eval(x))/y*100);
    grxfG4->SetPoint(n,x,(y-ffxG4->Eval(x))/y*100);
    ++n;
  }
  grxfG4->SetLineColor(kBlue);
  grxfN->SetLineColor(kGreen+2);
  grxfG4->SetMarkerColor(kBlue);
  grxfN->SetMarkerColor(kGreen+2);
  grxfG4->SetMarkerStyle(23);
  grxfN->SetMarkerStyle(8);
  mgxf->SetTitle("Felicia Residual vs. Current (UK Data - Parametrization)");
  mgxf->GetXaxis()->SetTitle("Current/300 A");
  mgxf->GetYaxis()->SetTitle("Diff: UK Data - Pol5 (%)");
  gPad->SetGrid();
  mgxf->Add(grxfG4);
  mgxf->Add(grxfN);
  mgxf->Draw("acp");
  tl2->Draw();
  cxf->SaveAs("Felicia_G4_UK_residual.pdf");
}
