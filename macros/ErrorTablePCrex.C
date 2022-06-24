{
  std::vector<std::string>name;
  name.push_back("$\\langle A_{zz}\\rangle$");
  name.push_back("Beam Trajectory");
  name.push_back("Foil Polarization");
  name.push_back("Dead Time");
  name.push_back("Charge Normalization");
  name.push_back("Leakage Currents");
  name.push_back("Laser Polarization");
  name.push_back("Accidentals");
  name.push_back("Current Dependence");
  name.push_back("Aperture Transmission");
  name.push_back("Null Asymmetry");
  name.push_back("July Extrapolation");
  
  std::map<std::string,double>prex_err, crex_err;
  prex_err["$\\langle A_{zz}\\rangle$"] = 0.202;
  prex_err["Beam Trajectory"] = 0.30;
  prex_err["Foil Polarization"]=0.628;
  prex_err["Dead Time"]=0.051;
  prex_err["Charge Normalization"]=0.0;
  prex_err["Leakage Currents"]=0.000;
  prex_err["Laser Polarization"]=0.100;
  prex_err["Accidentals"]=0.024;
  prex_err["Current Dependence"]=0.42;
  prex_err["Aperture Transmission"]=0.100;
  prex_err["Null Asymmetry"]=0.120;
  prex_err["July Extrapolation"]=0.235;

  crex_err["$\\langle A_{zz}\\rangle$"] = 0.156;
  crex_err["Beam Trajectory"] = 0.00;
  crex_err["Foil Polarization"]=0.571;
  crex_err["Dead Time"]=0.147;
  crex_err["Charge Normalization"]=0.01;
  crex_err["Leakage Currents"]=0.18;
  crex_err["Laser Polarization"]=0.06;
  crex_err["Accidentals"]=0.041;
  crex_err["Current Dependence"]=0.5;
  crex_err["Aperture Transmission"]=0.1;
  crex_err["Null Asymmetry"]=0.22;
  printf("\\begin{tabular}[lcc]\\hline\n");
  printf("Uncertainty & PREX-2 & CREX \\\\\\hline\n");
  double p_err=0, c_err=0;
  for (std::vector<std::string>::iterator it = name.begin();it!=name.end();++it){
    printf("%s & %0.2f & %0.2f \\\\\n",(*it).data(),prex_err[*it],crex_err[*it]);
    p_err += pow(prex_err[*it],2);
    c_err += pow(crex_err[*it],2);
  }
  printf("\\hline\n");
  printf("Total&%0.2f&%0.2f\\\\\\hline\n",sqrt(p_err),sqrt(c_err));
  printf("\\end{tabular}\n");

}
