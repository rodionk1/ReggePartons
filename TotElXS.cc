#include <fstream>
#include <iostream>
//#include <string> 
#include "stdlib.h"                                                                                                                                        
#include <iomanip>
#include "math.h"

int main(int argc, char *argv[])
{
  /* Builtin functions */
  double sin(double);
  const static int nc =50;
  const static int rap = 14;
  const static double rapstep = 1.0;
  const static double rapstart = 5.;
  /* Local variables */
  static double rint[nc], f[nc], h__;
  static double ampl[rap][nc], ampl2[rap][nc], sigmatot[rap][nc], ampl1n[rap][nc], sigma1n[rap][nc],sigmael[rap][nc]; 
  static double amplppbar[rap][nc], ampl2ppbar[rap][nc],sigmatotppbar[rap][nc],sigmaelppbar[rap][nc];
  static int i__;
  extern int qtsvd_c(double *, double *, int *, double *);
  static double bb[nc], bb1n[nc], x1;

  for (int i =0; i<nc; i++){
    for (int j=0;j<rap;j++){ampl[j][i]=0;ampl2[j][i]=0;};
  };

  /* testing on quadratic polynomials */
  /*
  static double testf[4];
  static double testx[4];
  static double testres[4];

  testx[0]=0.;
  testx[1]=0.5;
  testx[2]=0.75;
  testx[3]=1.;
  testf[0]=0.;
  testf[1]=0.25+0.5;
  testf[2]=0.5625+0.75;
  testf[3]=1.+1.;
  
  int testn=4;
  qtsvd_c(testres, testx, &testn, testf);
  std::cout << testres[0] << " " << testres[1] <<" " << testres[2]<< " "<<testres[3] << std::endl;
  */

  if(argc == 1){std::cout << "Please specify the input file."<<std::endl; return 0; }
  //  if(argc == 2){std::cout << "Missing input or output file. Please specify."<<std::endl; return 0; }

  int ibb;

  /* input file, amplitude p-p */
  for(int i=1; i<argc; i++){
  std::ifstream if1;
  if1.open(argv[i],std::ios::in);
  if (! if1)
    {
      std::cout << "The particle-particle amplitude file not found." << std::endl; return 0;
    };

  /* input file, amplitude Pomeron-p */
  /*  std::ifstream if2;
  if2.open(argv[2],std::ios::in);
  if (! if2)
    {
      std::cout << "The Pomeron-particle amplitude file not found." << std::endl; return 0;
    };
  
  */

  /* Reading data. */

  char s1[1000];
  char s;
  /*ampl*/
  if1.getline(s1,1000);
  double a;
  ibb=0;
  while(if1.get(s)){
    if1.putback(s);
    if1 >> bb[ibb]; 
    //    std::cout << s <<" " << bb[ibb] << std::endl;
    for (int irap=0; irap<rap; irap++){
      double a;
      if1 >> a;
      ampl[irap][ibb]+=a;    
    };
    ibb++;
     if1.getline(s1,1000);
  };
  /*a1n*/
  /************
  if2.getline(s1,100);
  int ibb1n=0;
  while(if2.get(s)){
    if2.putback(s);
    if2 >> bb1n[ibb1n];
    for (int irap=0; irap<10; irap++){
      if2 >> ampl1n[irap][ibb1n];};
    ibb1n++;
    if2.getline(s1,100);
  };
  *********/
  if1.close();
  };

  /*

  */


  for (int irap=0; irap<rap; irap++){
    for(int jb=0; jb<ibb; jb++){
      //   ampl[irap][jb]/=4.;
      ampl[irap][jb]/=(double) (argc-1);

    };
  };



  double pi = 3.14159265;


  double Deltaplus = -0.38;   // PDG -0.34 Volkovitsky -0.55 T-M -0.27 Poghosyan -0.3
  double alplus = 0.70*0.197*0.197; // Volkovitsky 0.70   T-M 0.70  Poghosyan 0.8  
  double betaplus2=6.01/(4*pi); // PDG 6.01                                                    
  double rplus = 3.0*0.197*0.197; // Volkovitsky 1.55 T-M 3.00 Poghosyan 0.918  

  double Deltaminus = -0.55;   // PDG -0.55 Volkovitsky -0.57 T-M -0.54 Poghosyan -0.6
  double alminus = 1.0*0.197*0.197; // Volkovitsky 1.0 T-M 1.0  Poghosyan 0.9 
  double betaminus2= 3.28/(4*pi); //PDG 3.28                                                                                                        
  double rminus = 10.0*0.197*0.197; // Volkovitsky 5.19 T-M 10.00 Poghosyan 0.945


  /* creating a^2 */
  for(int irap=0; irap<rap; irap++){
    for(int ii=0; ii<ibb; ii++){
      //      std::cout << ii<<" "<<  bb[ii] << std::endl;
      double y = rapstart+rapstep*irap;

      double regplus = betaplus2*exp(Deltaplus*y)/(2*alplus*y+2*rplus)*exp(-bb[ii]*bb[ii]/(4*alplus*y+4*rplus));
      double regminus = betaminus2*exp(Deltaminus*y)/(2*alminus*y+2*rminus)*exp(-bb[ii]*bb[ii]/(4*alminus*y+4*rminus));
      double amplP = ampl[irap][ii];

      ampl2[irap][ii]=((amplP+1)*(1 + regplus-regminus)-1)*((amplP+1)*(1+regplus-regminus)-1);
      ampl2ppbar[irap][ii]=((amplP +1)*(1+regplus+regminus)-1)*((amplP +1)*(1+regplus+regminus)-1);

      double repp = (-(cos(pi*(Deltaplus+1))-1)*regplus/sin(pi*(Deltaplus+1))+(cos(pi*(Deltaminus+1))-1)*regminus/sin(pi*(Deltaminus+1)));
      double reppbar = (-(cos(pi*(Deltaplus+1))-1)*regplus/sin(pi*(Deltaplus+1))-(cos(pi*(Deltaminus+1))-1)*regminus/sin(pi*(Deltaminus+1)));
      /* replus -- real part of the pp amplitude  */
      /* reminus -- real part of the ppbar amplitude */

      ampl2[irap][ii]+= repp*repp*(amplP+1)*(amplP+1);
      ampl2ppbar[irap][ii]+= reppbar*reppbar*(amplP+1)*(amplP+1);


      amplppbar[irap][ii]=(amplP + 1)*(1+regplus+regminus)-1;
      ampl[irap][ii] = (amplP + 1)*(1+regplus-regminus)-1;

      //      std::cout<< amplP << std::endl;


    };
  };


  /* eikonal amplitudes */
  // Ter-Martirosyan settings
  /*  double RR=0.36;
  double Delta = 0.12;
  double gamma = 0.08305126;
  double alp = 0.008538;
  */

  // Shabelski settings
  double RR=0.3513;
  double Delta = 0.139;
  double gamma = 0.06869;
  double alp = 0.00814989;
 

 std::cout << "#y  stot_pp  stot_ppbar  sel_pp  sel_ppbar" << std::endl;


  /* computing total cross section for all rapidities*/
  for (int irap = 0; irap<rap; irap++){
    static double bampl[nc], bamplppbar[nc];
    for (int ii=0; ii<ibb; ii++){
      bampl[ii]=4.*3.14159265*ampl[irap][ii]*bb[ii];
      bamplppbar[ii]=4.*3.14159265*amplppbar[irap][ii]*bb[ii];
      //      std::cout << bb[ii] <<" "<<bampl[ii] << std::endl;
    };
    double resint[nc],resinteik[nc];
qtsvd_c(resint, bb, &ibb, bampl);
 for (int ii=0; ii<nc; ii++){
   sigmatot[irap][ii] = resint[ii];  
     };
 qtsvd_c(resinteik, bb, &ibb, bamplppbar);
 for (int ii=0; ii<nc; ii++){
   sigmatotppbar[irap][ii] = resinteik[ii];
 };

   
 static double bampl2[nc],bampl2ppbar[nc];
  for (int ii=0; ii < ibb; ii++){
    bampl2[ii]=2.*3.14159265*ampl2[irap][ii]*bb[ii];
    bampl2ppbar[ii]=2.*3.14159265*ampl2ppbar[irap][ii]*bb[ii];
  };
  double resel[nc],reseleik[nc];
  qtsvd_c(resel, bb, &ibb, bampl2);
  qtsvd_c(reseleik, bb, &ibb, bampl2ppbar);
 for (int ii=0; ii<nc; ii++){
   sigmael[irap][ii] = resel[ii];  
   sigmaelppbar[irap][ii] = reseleik[ii];  
     };
 /***********
 static double bampl1n[nc];
 for (int ii=0; ii<ibb1n; ii++){
   bampl1n[ii]=ampl1n[irap][ii]*bb1n[ii];
 };
 double res1n[nc];
  qtsvd_c(res1n, bb, &ibb1n, bampl1n);
 
  for (int ii=0; ii<nc; ii++){
    sigma1n[irap][ii]=res1n[ii];    
  };
 ****************/

 std::cout << exp((rapstart+rapstep*irap)/2.) <<" " <<sigmatot[irap][ibb-1]*10 <<" " <<sigmatotppbar[irap][ibb-1]*10 <<" ";
    std::cout<<" " << sigmael[irap][ibb-1]*10 <<" " << sigmaelppbar [irap][ibb-1]*10 <<" ";
  //  std::cout <<" "<< pow(sigma1n[irap][ibb1n-1]/0.03,2)/(sigmatot[irap][ibb-1]-1.5*sigmael[irap][ibb-1]);
  //  std::cout <<" "<< exp(0.12*(1+2*irap))/(sigmatoteik[irap][ibb-1]-1.5*sigmaeleik[irap][ibb-1]) ;
  std::cout << std::endl;  
  }; std::cout << "#*" <<std::endl;



  return 0;
} /* main */
