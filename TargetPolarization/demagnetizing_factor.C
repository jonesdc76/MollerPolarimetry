{
  //Uses equations 1a and 1b from "Effective Demagnetizing Factors of 
  //Complicated Particle Mixtures", 2007, Ralph Skompski et al.
  gStyle->SetTitleW(0.9);
  TF1 *f = new TF1("f1","(x<0.999 ? pow(1/(1-x*x+0.000001)*(1-x/sqrt(1-x*x+0.000001)*acos(x)),-1):(x>1.001 ? pow(1/(x*x-1)*(x/sqrt(x*x-1)*acosh(x)-1),-1) : 1+2*x))",0.00001,15);
  f->SetLineColor(kBlue);
  f->SetLineWidth(3);
  f->SetNpx(500);
  TCanvas *cd = new TCanvas("cd","cd", 0, 0, 660, 500);
  cd->SetGrid();
  f->Draw();
  f->SetTitle("Ellipsoid of Rotation Demagnetizing Factor vs Axis Ratio");
  f->GetXaxis()->SetTitle("Ellipsoid Axis Ratio (R_{z}/R_{x})");
  f->GetYaxis()->SetTitle("Demagnetizing Factor");
  //cd->SaveAs("demagnetizing_factor.pdf");
}
