{
  double in4_2mm = 88.78, in4_2mm_err=0.27;
  double out4_neg2mm = 90.28, out4_neg2mm_err=0.29;
  double out4_2mm = 89.89, out4_2mm_err=0.27;
  double out4_0mm = (89.03+89.34)/2.0, out4_0mm_err=0.28/sqrt(2);//0.28
  double in4_0mm = 87.72, in4_0mm_err = 0.27;	       
  double rat1=in4_2mm/in4_0mm, rat2=0.5*(out4_2mm+out4_neg2mm)/out4_0mm;
  double err1 = sqrt(pow(in4_2mm_err/in4_2mm,2)+pow(in4_0mm_err/in4_0mm,2)),err2=sqrt(pow(out4_neg2mm_err/out4_neg2mm,2)/2.0+pow(out4_0mm_err/out4_0mm,2))/rat2;
  cout<<"Ratio IHWP IN  +/-2 mm  to center: "<<rat1<<"+/-"<<err1<<endl;
  cout<<"Ratio IHWP OUT +/-2 mm  to center: "<<rat2<<"+/-"<<err2<<endl;
  cout<<"Avergage Ratio +/-2 mm  to center: "<<(rat1/pow(err1,2)+rat2/pow(err2,2))/(pow(err1,-2)+pow(err2,-2))<<"+/-"<<1/sqrt(pow(err1,-2)+pow(err2,-2))<<endl;

}
