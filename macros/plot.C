{
  gStyle->SetOptFit(1111);
  TF1 *f = new TF1("f","pol0",0,1);
  TF1 *fh = new TF1("fh","pol0",0,1);
  TGraphErrors *grl = new TGraphErrors("lowfield.txt");
  grl->SetMarkerStyle(8);
  TGraphErrors *grh = new TGraphErrors("highfield.txt");
  grh->SetMarkerStyle(8);
  grh->SetMarkerColor(kRed);
  TGraphErrors *grlin = new TGraphErrors("lowfieldNeg.txt");
  int len = grl->GetN();
  for(int i=0;i<grlin->GetN();++i){
    double x, y;
    grlin->GetPoint(i,x,y);
    grl->SetPoint(len+i,x,y);
    grl->SetPointError(len+i, 0, grlin->GetErrorY(i));
    cout<<len+i<<" "<<x<<" "<<y<<" "<<grlin->GetErrorY(i)<<endl;
  }
  grlin->SetMarkerStyle(8);
  grlin->SetMarkerColor(kMagenta);
  TGraphErrors *grhin = new TGraphErrors("highfieldNeg.txt");
  len = grh->GetN();
  for(int i=0;i<grhin->GetN();++i){
    double x, y;
    grhin->GetPoint(i,x,y);
    grh->SetPoint(i+len,x,y);
    grh->SetPointError(i+len, 0, grhin->GetErrorY(i));
  }
  grhin->SetMarkerStyle(8);
  grhin->SetMarkerColor(kBlue);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(grl);
  mg->Add(grh);
  mg->Add(grlin);
  mg->Add(grhin);
  f->SetLineColor(kBlack);
  f->SetRange(17920,17990);
  mg->Draw("ap");
  grl->Fit(f);
  grlin->SetLineColor(kMagenta);
  cout<<"3.5 T polarization: "<<f->GetParameter(0)<<"+/-"<<f->GetParError(0)<<endl;
  grh->SetLineColor(kRed);
  fh->SetLineColor(kRed);
  grh->Fit(fh);
  mg->Draw("ap");
  grhin->SetLineColor(kBlue);
  //f->SetLineColor(kBlue);
  //grhin->Fit(f);
  cout<<"4.2 T polarization: "<<fh->GetParameter(0)<<"+/-"<<fh->GetParError(0)<<endl;
  mg->SetTitle("Measured Polarization for 3.5 and 4.0 T Target Coil Fields");
  mg->GetXaxis()->SetTitle("Run");
  mg->GetYaxis()->SetTitle("Polarization");
  double x[2] = {3.5,4.0},xe[2]={0,1};
  double y[2] = {88.3518, 89.0172}, ye[2]={0.1896,0.127};
  TGraphErrors *gr = new TGraphErrors(2,x,xe,y,ye);
}
