{
  TRandom r(0);
  TH1D *h = new TH1D("h","h",100,-0.3,0);
  TH1D *hpout = new TH1D("hpout","hpout",100,0.96,1.0);
  TH1D *hpin = new TH1D("hpin","hpin",100,0.96,1.0);
  TF1 *f;
  double eps = 0.000001, step = 0.01;
  for(int i=0;i<10000;++i){
    step = 0.01;
    double pout = r.Gaus(0.164,0.02), pin = r.Gaus(-0.053,0.02);
    TF1 f("f",Form("sqrt(1-pow(%f+x,2))-sqrt(1-pow(%f+x,2))-%f",pout,pin,r.Gaus(0.012,0.0014)),0,1);
    double dvw = 0, sign = (f.Eval(dvw)>0 ? -1.0:1.0), val = f.Eval(dvw), prev_val = val;

    while(abs(val)>eps){
      if(abs(val) > abs(prev_val)){
	step *= 0.5;
	sign *= -1.0;
      }
      dvw += step * sign;
      prev_val = val;
      val = f.Eval(dvw);
      //intf("%f, %f, %f\n",dvw, step, val);
    }
    if(i%1000==0)printf("%i  %f\n",i, dvw);
    hpout->Fill(sqrt(1-pow(pout+dvw,2)));
    hpin->Fill(sqrt(1-pow(pin+dvw,2)));
    h->Fill(dvw);
  }

  h->Draw();
new TCanvas;
hpin->SetLineColor(kRed);
hpout->Draw();
hpin->Draw("sames");

}
