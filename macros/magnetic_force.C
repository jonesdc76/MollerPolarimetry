{
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.05);
  double mu_Fe = 220;//emu/g
  double rho_Fe = 7.874;//g/cm^3
  double foil_thickness = 0.1;//in cm
  double emu_to_J_per_T = 1.0e-3;
  TCanvas *c = new TCanvas("c","c",0,0,1000,800);
  c->Divide(1,2);
  c->cd(1)->SetGrid();
  TGraph *gr = new TGraph("magnetic_field.dat","%lg%lg");
  for(int i=0;i<gr->GetN();++i){
    double x, y;
    gr->GetPoint(i,x,y);
    gr->SetPoint(i,x,y/1.0e4);
  }
  gr->SetTitle("Magnetic Field Along Beam Line Measured from Coil Center");
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(20);
  gr->SetLineColor(kBlue);
  gr->SetLineWidth(2);
  gr->Draw("acp");
  gr->GetXaxis()->SetTitle("Distance from Center Along Beamline (cm)");
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitle("B-field (T)");
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleOffset(0.5);
  c->cd(2)->SetGrid();
  TGraph *grf = new TGraph("magnetic_field.dat","%lg%*lg%lg");
  for(int i=0;i<grf->GetN();++i){
    double x, y;
    grf->GetPoint(i,x,y);
    x -= 0.635;
    y *= TMath::Pi()*foil_thickness*pow(2.54/2.0,2)*rho_Fe*mu_Fe*emu_to_J_per_T;
    grf->SetPoint(i,x,y);
  }
  grf->SetTitle(Form("Force on Foil vs Z along Beam (5 T field, 25.4 mm diameter %0.2f mm foil)", foil_thickness*10.0));
  grf->SetMarkerColor(kRed);
  grf->SetMarkerStyle(20);
  grf->SetLineColor(kRed);
  grf->SetLineWidth(2);
  grf->Draw("acp");
  grf->GetXaxis()->SetTitle("Distance from Center Along Beamline (cm)");
  grf->GetXaxis()->SetTitleSize(0.05);
  grf->GetYaxis()->SetTitle("Force (N)");
  grf->GetYaxis()->SetTitleOffset(0.5);
  grf->GetYaxis()->SetTitleSize(0.05);
  grf->SetLineColor(kRed);
  c->SaveAs("MagneticForceOnFoil.pdf");
}







