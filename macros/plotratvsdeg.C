{

   auto c1 = new TCanvas("c1","c1",800, 500);
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Nomalized Ratio 1F/2F vs B-Field");
  TGraph *gr[6];
  int col[6] = {kBlack, kBlue, kGreen+2, kRed, kMagenta, kOrange};
  int shape[6] = {20,21,22,23,33,34};
  TString name[6] = {"+2 deg","+1 deg","+0.5 deg","-0.5 deg","-1 deg","-2 deg"};
  TString draw = "%lg";
  TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
  for(int i=0;i<6;++i){
    TString dr = draw + "%lg";
    gr[i] = new TGraph("polvsdegnorm.dat",dr.Data());
    gr[i]->SetName(name[i].Data());
    draw += "%*lg";
    gr[i]->SetMarkerColor(col[i]);
    gr[i]->SetMarkerStyle(shape[i]);
    gr[i]->SetMarkerSize(0.8);
    mg->Add(gr[i]);
    leg->AddEntry(gr[i],name[i].Data());
  }

  mg->Draw("alp");
  mg->GetXaxis()->SetTitle("B-field (T)");
  mg->GetYaxis()->SetTitle("1F/2F Normalized to 4T");
  leg->Draw();

  new TCanvas;
  double x[6] = {2,1,0.5,-0.5,-1,-2};
  double y[6] = {0.8211,0.8372,0.8488,0.8438,0.8360,0.8179};
  TGraph *gr1 = new TGraph(6,x,y);
  gr1->SetMarkerStyle(8);
  gr1->Draw("ap");
  gr1->GetXaxis()->SetTitle("Foil Readback Angle (#circ)");
  gr1->GetYaxis()->SetTitle("Normalized 1F/2F @ 1.82 T");
  gPad->Update();
}
