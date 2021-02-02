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
  double cos(double);
  const static int nc =50;
  const static int rap = 6;
  double rapstep = 1.0;
  double rapstart = 11.;
  /* Local variables */
  static double rint[nc], f[nc], h__;
  static double amplppbar[rap][nc], ampl2ppbar[rap][nc], sigmatotppbar[rap][nc], ampl1n[rap][nc], sigma1n[rap][nc],sigmaelppbar[rap][nc]; 
  static double ampl[rap][nc], ampl2[rap][nc],sigmatot[rap][nc],sigmael[rap][nc];
  static double chi11[rap][nc], chi12[rap][nc], chi22[rap][nc];
  static double amplsd[rap][nc],sigmasd[rap][nc];
  static double amplsd2[rap][nc],sigmasd2[rap][nc];
  static int i__;
  extern int qtsvd_c(double *, double *, int *, double *);
  static double bb[nc], bb1n[nc], x1;


  /* input file, amplitude p-p */
  std::ifstream if1;
  if1.open(argv[1],std::ios::in);
  if (! if1)
    {
  std::cout <<"Usage:  ./executable  amplitude.dat SDampl.dat" << std::endl ;
      std::cout << "The particle-particle amplitude file not found." << std::endl; return 0;
    };
 


  /* Reading data. */

  char s1[1000];
  char s;
  /*ampl*/
  if1.getline(s1,1000);
  double as;
  int ibb=0;
  while(if1.get(s)){
    if1.putback(s);
    if1 >> bb[ibb];
    //    std::cout << s <<" " << bb[ibb] << std::endl;                                                                                                      
    for (int irap=0; irap<rap; irap++){
      if1 >> as; 
      ampl[irap][ibb]=as;};
    ibb++;
    if1.getline(s1,1000);
  };



  /* input file, amplitude SD */
 
  ibb=0;

  std::ifstream if2;
  if2.open(argv[2],std::ios::in);
  if (! if2)
    {
  std::cout <<"Usage:  ./executable  amplitude.dat SDampl.dat" << std::endl ;
      std::cout << "The SD amplitude file #1 is not found." << std::endl; return 0;
    };


  /* Reading data. */

  char s1SD[1000];
  char sSD;
  /*ampl*/
  if2.getline(s1SD,1000);

  while(if2.get(sSD)){
    if2.putback(sSD);
    //    if1 >> bb[ibb];
    if2 >> as;
    //    std::cout << as <<" " << bb[ibb] << std::endl;                                                                     
    for (int irap=0; irap<rap; irap++){
      if2 >> as;
      amplsd[irap][ibb]=0; 
    amplsd[irap][ibb]+=as;};
    ibb++;
    if2.getline(s1SD,1000);
  };





  /* input file, amplitude SD */
  /*
  */


  /* input file, amplitude SD */
  /*

  */



  // Shabelski settings
  double RR=0.3513;
  double Delta = 0.145;
  double gamma = 0.06869;
  double alp = 0.007;

  double pi = 3.14159265;
 
  double R1 = 0.335;
  double R2 = 0.335;
  double a = 0.018;
  double eps = 3.14159265*a*a;
  double Nbar = 34.5;
  double eta = 0.75;
  double N1 = Nbar*(1+eta);
  double N2 = Nbar*(1-eta);

   
  double Deltaplus = -0.38; // PDG -0.34 Volkovitsky -0.55 T-M -0.27 Poghosyan -0.3
  double alplus = 0.70*0.197*0.197; // Volkovitsky 0.70   T-M 0.70  Poghosyan 0.8
  double betaplus2=6.01/(4*pi); // PDG 6.01
  //  double betaplus2=0.; //15*0.197*0.197;
  double rplus = 3.0*0.197*0.197; // Volkovitsky 1.55 T-M 3.00 Poghosyan 0.918

  double Deltaminus = -0.55;  // PDG -0.55 Volkovitsky -0.57 T-M -0.54 Poghosyan -0.6
  double alminus = 1.0*0.197*0.197; // Volkovitsky 1.0 T-M 1.0  Poghosyan 0.9 
  double betaminus2= 3.28/(4*pi); //PDG 3.28
  //  double betaminus2= 0.; //1.8*0.197*0.197;
  double rminus = 10*0.197*0.197; // Volkovitsky 5.19 T-M 10.00 Poghosyan 0.945


  for (int irap =0; irap<rap; irap++){
    for(int ii=0;ii<ibb; ii++){
      //      std::cout <<      bb[ii] << std::endl;
      double y = rapstart+irap*rapstep;

      double regplus = betaplus2*exp(Deltaplus*y)/(2*alplus*y+2*rplus)*exp(-bb[ii]*bb[ii]/(4*alplus*y+4*rplus));
      double regminus = betaminus2*exp(Deltaminus*y)/(2*alminus*y+2*rminus)*exp(-bb[ii]*bb[ii]/(4*alminus*y+4*rminus));
      double amplP = ampl[irap][ii]; //  ((1-exp(-chi11[irap][ii]))+2*(1-exp(-chi12[irap][ii]))+(1-exp(-chi22[irap][ii])))/4.   ;
 
     //      ampl[irap][ii] = amplP;
      //      std::cout<<bb[ii]<<" " << ampleik[irap][ii] << " " << std::endl;
      ampl2[irap][ii]=((amplP+1)*(1 + regplus-regminus)-1)*((amplP+1)*(1+regplus-regminus)-1);
      ampl2ppbar[irap][ii]=((amplP +1)*(1+regplus+regminus)-1)*((amplP +1)*(1+regplus+regminus)-1);
      
double repp = (-(cos(pi*(Deltaplus+1))-1)*regplus/sin(pi*(Deltaplus+1))+(cos(pi*(Deltaminus+1))-1)*regminus/sin(pi*(Deltaminus+1)));
double reppbar = (-(cos(pi*(Deltaplus+1))-1)*regplus/sin(pi*(Deltaplus+1))-(cos(pi*(Deltaminus+1))-1)*regminus/sin(pi*(Deltaminus+1)));
     

/* replus -- real part of the pp reggeon exchange amplitude  */
/* reminus -- real part of the ppbar reggeon exchange amplitude */

      ampl2[irap][ii]+= repp*repp*(amplP+1)*(amplP+1); 
      ampl2ppbar[irap][ii]+= reppbar*reppbar*(amplP+1)*(amplP+1);
      amplppbar[irap][ii]=(amplP + 1)*(1+regplus+regminus)-1;
      ampl[irap][ii] = (amplP + 1)*(1+regplus-regminus)-1;

  //      std::cout << 2*((1+2*(regplus-regminus)+repp*repp + (regplus-regminus)*(regplus-regminus))*amplsd[irap][ii] + 
      //		      (repp*repp + (regplus-regminus)*(regplus-regminus))*(1+2*amplP) +2*amplP*(regplus-regminus)   -ampl2[irap][ii]) << " " ;


      
 amplsd2[irap][ii]=2*(amplsd[irap][ii]-amplP*amplP)*(1+repp*repp + (regplus-regminus)*(regplus-regminus)+ 2*(regplus-regminus));


 // std::cout << amplsd2[irap][ii] << " " << amplsd[irap][ii]<<" "<< amplP*amplP <<" " << (1+repp*repp + (regplus-regminus)*(regplus-regminus)+ 2*(regplus-regminus)) << " *  ";

    };   std::cout << std::endl;
  };

  std::cout<< "# + " << Deltaplus << " " << rplus << " "<< alplus << " " << betaplus2 << std::endl;
  std::cout<< "# - " << Deltaminus << " " << rminus << " "<< alminus << " " << betaminus2 << std::endl;
//  std::cout<< "# P " << Delta << " " << R1 <<" "<< R2 << " "<< alp << " " << N1 <<" " << N2 << std::endl; 
 std::cout << "#y  s_tot  s_el" << std::endl;


  /* computing total cross section for all rapidities*/
  for (int irap = 0; irap<rap; irap++){
    static double bampl[nc], bampleik[nc];
    for (int ii=0; ii<ibb; ii++){
      bampl[ii]=4.*3.14159265*amplppbar[irap][ii]*bb[ii];
      bampleik[ii]=4.*3.14159265*ampl[irap][ii]*bb[ii];
      //      std::cout << bb[ii] <<" "<<bampl[ii] << std::endl;
    };
    double resint[nc],resinteik[nc];
qtsvd_c(resint, bb, &ibb, bampl);
 for (int ii=0; ii<nc; ii++){
   sigmatotppbar[irap][ii] = resint[ii];  
     };
 qtsvd_c(resinteik, bb, &ibb, bampleik);
 for (int ii=0; ii<nc; ii++){
   sigmatot[irap][ii] = resinteik[ii];
 };

   
 static double bampl2[nc],bampl2eik[nc];
  for (int ii=0; ii < ibb; ii++){
    bampl2[ii]=2.*3.14159265*ampl2ppbar[irap][ii]*bb[ii];
    bampl2eik[ii]=2.*3.14159265*ampl2[irap][ii]*bb[ii];
  };
  double resel[nc],reseleik[nc];
  qtsvd_c(resel, bb, &ibb, bampl2);
  qtsvd_c(reseleik, bb, &ibb, bampl2eik);
 for (int ii=0; ii<nc; ii++){
   sigmaelppbar[irap][ii] = resel[ii];  
   sigmael[irap][ii] = reseleik[ii];  
     };
 /*
 static double bamplsd[nc];
 for (int ii=0; ii < ibb; ii++){
   bamplsd[ii]=2.*3.14159265*amplsd[irap][ii]*bb[ii];
 };
 double reselsd[nc];
 qtsvd_c(reselsd, bb, &ibb, bamplsd);
 for (int ii=0; ii<nc; ii++){
   sigmasd[irap][ii] = reselsd[ii];
 };
 */

 static double bamplsd2[nc];
 for (int ii=0; ii < ibb; ii++){
   bamplsd2[ii]=2.*3.14159265*amplsd2[irap][ii]*bb[ii];
 };
 double reselsd2[nc];
 qtsvd_c(reselsd2, bb, &ibb, bamplsd2);
 for (int ii=0; ii<nc; ii++){
   sigmasd2[irap][ii] = reselsd2[ii];
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

 std::cout << exp((rapstart +rapstep*irap)/2.) <<" " <<sigmatot[irap][ibb-1]*10 <<" " <<sigmatotppbar[irap][ibb-1]*10 <<" ";
 std::cout<<" " << sigmael[irap][ibb-1]*10 <<" " << sigmaelppbar[irap][ibb-1]*10 <<" "<<sigmasd2[irap][ibb-1]*10;// << " "<<sigmasd2[irap][ibb-1]*10;;
  //  std::cout <<" "<< pow(sigma1n[irap][ibb1n-1]/0.03,2)/(sigmatot[irap][ibb-1]-1.5*sigmael[irap][ibb-1]);
  //  std::cout <<" "<< exp(0.12*(1+2*irap))/(sigmatoteik[irap][ibb-1]-1.5*sigmaeleik[irap][ibb-1]) ;
  std::cout << std::endl;  
  }; std::cout << "#*" <<std::endl;



  return 0;
} /* main */
