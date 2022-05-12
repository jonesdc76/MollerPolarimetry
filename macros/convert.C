#include<iostream>
#include<fstream>
#include <stdio.h>

using namespace std;
int convert(){
  ifstream f;
  FILE *fout = fopen("postcovid.csv","w");
  f.open("postcovid.txt");
  char a[100],b[100],c[100],d[100],e[100];
  f>>a>>b>>c>>d;
  fprintf(fout,"%s,%s,%s,%s\n", a,b,c,d);  
  while (!f.eof()){
    f>>a>>b>>c>>d>>e;
    fprintf(fout,"%s %s, %s, %s, %s\n", a,b,c,d,e);  
  }
  return 0;
}
