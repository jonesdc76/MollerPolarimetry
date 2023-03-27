{
  //Choose the tabletype and set the beam current before running.
  //Produces table of values and errors for target foil polarization
  
  bool tabletype = 0;//two table styles are available
  double current = 1.0;//current in muA
  const double Z = 26, density = 55.845, bohr = 9.274009994e-21, gs = 2.002319304, avogadro = 6.0221409e23, TRisePerMuA=11.5;
  double dMsdT = -0.024;//slope of temp correction emu/g/degC
  double dT = TRisePerMuA*current;
  double Ms294 = 218.044, Ms294e = 0.436;//emu/g @T=294, B=4T
  double Tcorr = dT*dMsdT;//temp correction in emu/g/degC
  double Tcorre = 0.3*Tcorr;//error in temp correction in emu/g/degC
  //cout<<Tcorre<<endl;
  double MsHot = Ms294 + Tcorr, MsHote = sqrt(Ms294e*Ms294e+Tcorre*Tcorre);
  if(tabletype){
    printf("\\begin{tabular}{|l|l|l|c|}\\hline\n");
    printf("Quantity&Value&Error&Unit\\\\\\hline\n");
    printf("Saturation magnetization $M_s$(T=294 K) &%0.2f&%0.2f&emu/g\\\\\n",Ms294,Ms294e);
  }else{
    printf("\\begin{tabular}{|l|l|l|c|}\\hline\n");
    printf("Quantity&T=294 K&T=%i K&Unit\\\\\\hline\n",(int)round(294+dT));
    printf("Saturation magnetization $M_s$ &%0.2f(%i)&%0.2f(%i)&emu/g\\\\\n",Ms294,(int)round(Ms294e*100),MsHot,(int)round(MsHote*100));
  }
  Ms294*=density/bohr/avogadro;//muB/atom @T=294, B=4T
  Ms294e*=density/bohr/avogadro;//muB/atom @T=294, B=4T
  MsHot*=density/bohr/avogadro;//muB/atom @T=294, B=4T
  MsHote*=density/bohr/avogadro;//muB/atom @T=294, B=4T
  if(tabletype)
    printf("Saturation magnetization $M_s$(T=313 K) &%0.4f&%0.4f&$\\mu_B$/atom\\\\\n",MsHot,MsHote);
  else
    printf("Saturation magnetization $M_s$&%0.4f(%i)&%0.4f(%i)&$\\mu_B$/atom\\\\\n",Ms294,(int)round(Ms294e*1e4),MsHot,(int)round(MsHote*1e4));
    
  double gp = 1.9206, gpe=0.0019;
  double forb =  (gs-gp)/gp/(gs-1);
  double forbe = 2/gp/gp*gpe;
  double Mspin294 = Ms294*(1-forb);
  double Mspin294e = sqrt(pow((1-forb)*Ms294e,2)+pow(Ms294*forbe,2));
  double MspinHot = MsHot*(1-forb);
  double MspinHote = sqrt(pow((1-forb)*MsHote,2)+pow(MsHot*forbe,2));
  if(tabletype){
    printf("$g^{\\prime}$&%0.4f&%0.4f&$-$\\\\\n",gp,gpe);
    printf("Orbital fraction: $\\frac{M_{L}}{M_{\\rm tot}}=\\frac{g_{S}-g^\\prime}{g^\\prime(g_{S}-1)}$&%0.4f&%0.4f&$-$\\\\\n",forb,forbe);
    printf("Spin component: $M_S\\left(1-\\frac{M_{L}}{M_{\\rm tot}}\\right)$&%0.4f&%0.4f&$\\mu_B$/atom\\\\\n",MspinHot,MspinHote);
    printf("Avg. electron magnetization (T=313 K)&%0.6f&%0.6f&$\\mu_B$\\\\\n",MspinHot/Z,MspinHote/Z);
    printf("Avg. electron polarization (T=313 K)&%0.6f&%0.6f&$-$\\\\\\hline\n",MspinHot/Z/gs*2,MspinHote/Z/gs*2);
  }else{
    printf("$g^{\\prime}$&%0.4f(%i)&%0.4f(%i)&$-$\\\\\n",gp,(int)round(gpe*1e4),gp,(int)round(gpe*1e4));
    printf("Orbital fraction: $\\frac{M_{L}}{M_{\\rm tot}}=\\frac{g_{S}-g^\\prime}{g^\\prime(g_{S}-1)}$&%0.4f(%i)&%0.4f(%i)&$-$\\\\\n",forb,(int)round(forbe*1e4),forb,(int)round(forbe*1e4));
    printf("Spin component: $M_S\\left(1-\\frac{M_{L}}{M_{\\rm tot}}\\right)$&%0.4f(%i)&%0.4f(%i)&$\\mu_B$/atom\\\\\n",Mspin294,(int)round(Mspin294e*1e4),MspinHot,(int)round(MspinHote*1e4));
    printf("Average electron magnetization&%0.5f(%i)&%0.5f(%i)&$\\mu_B$\\\\\n",Mspin294/Z,(int)round(Mspin294e/Z*1e5),MspinHot/Z,(int)round(MspinHote/Z*1e5));
    printf("Average electron polarization&%0.6f(%i)&%0.5f(%i)&$-$\\\\\\hline\n",Mspin294/Z/gs*2,(int)round(Mspin294e/Z/gs*2*1e5),MspinHot/Z/gs*2,(int)round(MspinHote/Z/gs*2*1e5));
  }
  printf("\\end{tabular}\n");
}
