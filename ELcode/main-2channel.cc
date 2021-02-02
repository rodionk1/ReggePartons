/*******************
Computation of the amplitude with quasieikonal initial conditions.
*****************/
#include <fstream>
#include <iostream>
//#include <string>
#include <iomanip>
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "RGauss.hh"
#include "string.h"

#include "parton.hh"
#include "parton.cc"
#include "pstack.hh"
#include "pstack.cc"
#include "pbuffer.hh"
#include "pbuffer.cc"


int poisson(double );

double makeT(double *, int,int);
unsigned int graycode (unsigned int);
unsigned long int binom(unsigned int n, unsigned int k);
//double RandGauss();




int main(int argc, char *argv[])
{
  int Nevol = 20;
  int Nevent = 10000;
  double bb = 0.;
  double bmin = 0.;
  double bstep = 0.;
  int Nb = 1;
  int Niteration = 1;
  double Cquasi = 1.;
  const int Narray = 10;
  int CNlong = 100;
  double rapstep = 1.;
  double rapstart = 1.;
  double skscale = 1.;
  long double egapDelta = 0.;
  long double egapalp =0.;

  double lam,m1,m2,nu,alp;
  long double aa;
  long double Nbar1, R1, Pch1, Nbar2, R2, Pch2;
 long double egap=0.;
  double params[6];
  if(argc == 1){std::cout << "Please specify the input and output files."<<std::endl; return 0; }
  if(argc == 2){std::cout << "Missing input or output file. Please specify."<<std::endl; return 0; }

  /* input file */
  std::ifstream if1;
  if1.open(argv[1],std::ios::in);
  if (! if1)
    {
      std::cout << "The input file not found." << std::endl; return 0;
    };

  char s1[100];
  char s;

  while(if1.get(s)){
    if1.putback(s);
    while((if1.get(s))&&(!(s=='#'))){
      /*if not a comment, define what is it*/
      if1.putback(s);
      if1 >> s1;
      int num;
      if(strcmp(s1,"lam")==0) num=1;
      if(strcmp(s1,"m1")==0) num=2;
      if(strcmp(s1,"nu")==0) num=3;
      if(strcmp(s1,"m2")==0) num=4;
      if(strcmp(s1,"alp")==0) num=5;
      if(strcmp(s1,"R1")==0) num=6;
      if(strcmp(s1,"aa")==0) num=7;
      if(strcmp(s1,"Nbar1")==0) num=8;
      if(strcmp(s1,"Nevent")==0) num=9;
      if(strcmp(s1,"bb")==0) num=10;
      if(strcmp(s1,"Nevol")==0) num=11;
      if(strcmp(s1,"Nlong")==0) num =13;
      if(strcmp(s1,"rapstep")==0) num =14;
      if(strcmp(s1,"bmin")==0) num=15;
      if(strcmp(s1,"bstep")==0) num=16;
      if(strcmp(s1,"Nb")==0) num=17;
      if(strcmp(s1,"Niteration")==0) num=18;
      if(strcmp(s1,"sqrtkscale")==0) num=19;
      if(strcmp(s1,"rapstart")==0) num=20;
      if(strcmp(s1,"R2")==0) num =21;
      if(strcmp(s1,"Nbar2")==0) num =22;      
      if(strcmp(s1,"Pch1")==0) num =23;      
      if(strcmp(s1,"Pch2")==0) num =24;      
      if(strcmp(s1,"egap")==0) num =25;          
      if(strcmp(s1,"egapDelta")==0) num =26;          
      if(strcmp(s1,"egapalp")==0) num =27;          
      //      std::cout << num << std::endl;                                                                                                                 
      switch(num){
      case 1: if1 >> params[4]; lam = params[4]; break;
      case 2: if1 >> params[5]; m1 = params[5]; break;
      case 3: if1 >> params[1]; nu = params[1]; break;
      case 4: if1 >> params[0]; m2 = params[0]; break;
      case 5: if1 >> params[3]; alp = params[3]; break;
      case 6: if1 >> R1; break;
      case 7: if1 >> params[2]; aa = params[2]; break;
      case 8: if1 >> Nbar1; break;
      case 9: if1 >> Nevent; break;
      case 10: if1 >> bb; break;
      case 11: if1 >> Nevol; break;
      case 12: if1 >> Cquasi; break;
      case 13: if1 >> CNlong; break;
      case 14: if1 >> rapstep; break;
      case 15: if1 >> bmin; break;
      case 16: if1 >> bstep; break;
      case 17: if1 >> Nb; break;
      case 18: if1 >> Niteration; break;
      case 19: if1 >> skscale; break;
      case 20: if1 >> rapstart; break;
      case 21: if1 >> R2; break;
      case 22: if1 >> Nbar2; break;
      case 23: if1 >> Pch1; break;
      case 24: if1 >> Pch2; break;
      case 25: if1 >> egap; break;
      case 26: if1 >> egapDelta; break;
      case 27: if1 >> egapalp; break;
      }; if1.getline(s1,100);
    }; if1.getline(s1,100);
  };

  std::cout << "Parameters read off successfully." << std::endl;

  //  double Cpart=sqrt(Cquasi);
  //double NN=Nbar*Cpart; // "part" from parton -- influence of quasieikonal constant on the parton level          
  int Nmax= (int) 20*Nbar1;
  double kscale = skscale*skscale;
  long double pi = 3.14159265358979323846;
  //  double pi=3.14159265 ; 
 //  std::cout << NN << std::endl;

  srand(time(NULL)+1);
  /*  
  pstack * stack2 = new pstack();
  parton * p2 = new parton(0.,0.);
  stack2->put(*p2);
  delete p2;
  for (int i=1; i<=Nevol; i++){
    stack2->evolve(1);};
 

  */
  //  return 0; 
  
  double ampl[Nevol+1][Nb];
  int Nlong[Nevol+1][Nb];
  double ampl2[Nevol+1][Nb];
  int Nlong2[Nevol+1][Nb];
  for (int i=0; i<=Nevol; i++){
    for(int j=0; j<Nb;j++){
      ampl[i][j]=0; ampl2[i][j]=0; Nlong[i][j]=0; Nlong2[i][j]=0;
  };};
  
  pbuffer* buffer1; 

  /* //redundant  
  int Np2[Niteration];
  double xp2[Niteration][Nmax+1];
  double yp2[Niteration][Nmax+1];
  */

  /* now the same for cutting at the edge */


  for(int i=1; i<=Nevent;i++){
  long double eps =aa;
  //      std::cout << "eps=aa="<<aa<<std::endl;
  buffer1 = new pbuffer();
  pstack * stack1;
  parton * p1;  double xp,yp;
  // int n1 = poisson(Nbar);
  //  int  n1=1; // checking f11
  //  std::cout << n1 << " partons to put in the stack" <<std::endl;
  /*generating partons according to 2-channel eikonal*/
  int n1 =0;
  double  r=(1.*rand())/RAND_MAX; // defining the channel
  
  int nch;
  if (r<Pch1) nch=1; else nch=2;
  switch(nch){
  case 1:
    /*number of partons and positions for channel 1*/
     n1 = poisson(Nbar1);
    //    n1=1771;
   //  std::cout<< "channel1, " << n1 << " counts"<< std::endl;
  for (int ip=1; ip<=n1; ip++){
    xp = R1*RandGauss();
    yp = R1*RandGauss();
    p1 = new parton(xp,yp);
     stack1 = new pstack(&params[0]);
     stack1->put(*p1);
     //  std::cout << " adding parton to stack" <<std::endl;
  delete p1;
     //  std::cout << "parton is added into the stack, now removed"<<std::endl;
    buffer1->put(*stack1);
     //  std::cout << " putting stack into the buffer"<< std::endl;
  delete stack1;
  };  
    break;

  case 2:
    /*number of partons and their positions for chanel 2*/
       n1 = poisson(Nbar2);
       //  n1=1771;
   //  std::cout << "channel2, " << n1 <<" counts" <<std::endl;
  for (int ip=1; ip<=n1; ip++){
    xp = R2*RandGauss();
    yp = R2*RandGauss();
    p1 = new parton(xp,yp);
     stack1 = new pstack(&params[0]);
     stack1->put(*p1);
     //  std::cout << " adding parton to stack" <<std::endl;
  delete p1;
     //  std::cout << "parton is added into the stack, now removed"<<std::endl;
    buffer1->put(*stack1);
     //  std::cout << " putting stack into the buffer"<< std::endl;
  delete stack1;
  };
    break;
  };

                                                                                         
 
  /* n1 generated*/
   
  /* Instead of stack 2 which is not evolved we use a set of Gaussian profiles */

  //   return 0;

    if((i%10)==0){std::cout << i <<" out of "<< Nevent <<std::endl;};
    if(((i-1)%100)==0&&i>1){
  std::cout << "#b, y  ";
  for (int iev=1; iev<=Nevol; iev++){ std::cout << rapstart+(iev-1)*rapstep <<" ";}
   std::cout << std::endl; 
  for (int j=0; j<Nb; j++){
    std::cout << bmin+j*bstep <<" ";
    for (int iev =1; iev<=Nevol; iev++){
    //    of1 <<1+(iev-1)*rapstep <<" "<<ampl[iev] <<" "<<Nlong[iev]<<" "<< ampl2[iev]<<" "<<Nlong2[iev]<<std::endl;
      //    of1 <<1+(iev-1)*rapstep <<" "<< ampl2[iev]<<" "<<Nlong2[iev]<<std::endl;
      std::cout << ampl2[iev][j]/(i-1) << " "; 
    }; std::cout <<std::endl;
  };
 };
 //    std::cout << i << std::endl;

         buffer1->fullevolve(1.0*rapstart-egap);

       //       std::cout << "buffer 1 evolved" << std::endl;     


    for (int iev = 1;iev<=Nevol;iev++){

      //    std::cout << "   " << iev << " out of " << Nevol << std::endl;
    int count1 = 0; 
    int count2 = 0;    
     
      for  (int ic=1; ic<=buffer1->GetCount(); ic++){
	//	std::cout << ic <<std::endl;
	pstack* st1 = buffer1->read(ic);
	count1+=st1->GetCount();
      };

      
      //      std::cout << count1<<" and " << count2 <<" partons in the buffer, making T" << std::endl;

      const int ccount1 = count1;

      double astep;
      //      std::cout << "Evolution done, y = "<<rapstart+(iev-1)*rapstep <<"; "<< count1<< "partons"<< std::endl;
     
      /* full calculation */
      int flaglong=0;
      

      /* Calculation of the contribution to the amplitude, Beginning. */
   for (int iteration =0 ; iteration<Niteration; iteration++){
     /* Iterations for the buffer 2. */
     
     for(int ib=0; ib<Nb; ib++){
       double bcurrent = bmin+ ib*bstep;
       //       std::cout << iev << " " << bcurrent << std::endl ;     

      //      std::cout << "Building T matrix, count1=" << count1<<" count2="<<count2<<std::endl;

   if(ccount1==0){astep = 0; }else {
 long double G1[ccount1]; // CHANNEL 1
 long double G2[ccount1]; // CHANNEL 2
   int ic1=0;
   for (int iT = 1; iT<=buffer1->GetCount(); iT++){
     pstack* stack1 = buffer1->read(iT);
     for(int iTT=1; iTT<=stack1->GetCount(); iTT++){
	 p1 = stack1->read(iTT);
	 //	 std::cout << "Parton 1 read off" <<std::endl;
	 /*DO NOT DELETE STACK1 AND STACK2 -- THEY ARE _READ_ FROM BUFFER AND MUST STAY THERE*/
       double x = p1->getX(1) +bcurrent;
       double y = p1->getX(2);
       //    std::cout << x<<","<<y <<" ";
G1[ic1]= exp((long double)-(x*x+y*y)/(4*egapalp*egap+2*R1*R1));
G2[ic1]= exp((long double)-(x*x+y*y)/(4*egapalp*egap+2*R2*R2)); 

//G1[ic1]= 1./(2.*pi*(2*egapalp*egap+R1*R1))*exp((long double) egap*egapDelta)*exp((long double)-(x*x+y*y)/(4*egapalp*egap+2*R1*R1));
//G2[ic1]= 1./(2.*pi*(2*egapalp*egap+R2*R2))*exp((long double) egap*egapDelta)*exp((long double)-(x*x+y*y)/(4*egapalp*egap+2*R2*R2)); 
       //  std::cout <<T[ic1][ic2] << " ";        
       ic1++;
     };
   };
   
   int ccount = ic1; // number of partons in the stack
   //   std::cout << "G built" << "  ccount "<< ccount << " ccount1 "<< ccount1 <<std::endl;

   /* Now exclude zero lines and columns */
    /* Arrays of powers of G1 and G2*/   
   /* REDUNDANT
   double xk1[ccount][ccount+1]; // CHANNEL 1
   double xk2[ccount][ccount+1]; // CHANNEL 2   
   std::cout << "xxx" << std::endl;
   for (int ic=0; ic<ccount; ic++){
     xk1[ic][0]=1.0;      xk2[ic][0]=1.0; 
     std::cout << ic <<" " << G1[ic] << " " << G2[ic] << std::endl;
     for (int jc =1; jc<ccount+1; jc++){
       xk1[ic][jc]=xk1[ic][jc-1]*G1[ic];  
       xk2[ic][jc]=xk2[ic][jc-1]*G2[ic];
     };
     };*/
   //      std::cout << "a" <<std::endl; 
  long double sum1[ccount+1],sumk1[ccount+1]; // CHANNEL 1
  long double sum2[ccount+1],sumk2[ccount+1]; // CHANNEL 2
   sum1[0]=1; sumk1[0]=1;
   sum2[0]=1; sumk2[0]=1;
   for (int k=1; k<=ccount; k++){sum1[k]=0; sum2[k]=0;};

   for (int jc=1; jc<=ccount; jc++){ // 'jc' is power
     sumk1[jc]=0;     sumk2[jc]=0;
     for (int ic=0; ic<ccount; ic++){ // 'ic' is number 
       sumk1[jc]+=pow(G1[ic],jc);// xk1[ic][jc];       
       sumk2[jc]+=pow(G2[ic],jc);// xk2[ic][jc];      
     };
     //  std::cout<< sumk1[jc]<<" " << sumk2[jc];
     //  std::cout<< std::endl;
   };
   //  std::cout <<"b" << std::endl;
   for (int ic=0; ic<ccount; ic++){ // 'ic' is number                                                                                                      
     sum1[1]+=G1[ic];//    xk1[ic][1];  
     sum2[1]+=G2[ic];//    xk2[ic][1];
   };
   //  std::cout << "c" << std::endl;
  
   for (int k=1; k<=ccount-1; k++){
     for (int l=k; l>=0; l--){
       sum1[k+1]+=pow(-1,k-l)*sum1[l]*sumk1[k+1-l]; sum2[k+1]+=pow(-1,k-l)*sum2[l]*sumk2[k+1-l];
     };
     sum1[k+1]/=(k+1);     sum2[k+1]/=(k+1);
   };

   // std::cout << "d" << std::endl;
  
     long double fac1 = Nbar1*pi*eps*eps/(2.*pi*(2*egapalp*egap+R1*R1))*exp((long double) egap*egapDelta);
     long double fac2 = Nbar2*pi*eps*eps/(2.*pi*(2*egapalp*egap+R2*R2))*exp((long double) egap*egapDelta);
     //     std::cout << fac1 <<" " << fac2 << std::endl;
  long double result=0;
  /* int kmax=50;   
     if (kmax>ccount) kmax=ccount;
     for (int k=1; k<=kmax; k++){  */
     for (int k=1; k<=ccount; k++){
    long double a = pow(-1,k-1)*(Pch1*pow(fac1,k)*sum1[k] + Pch2*pow(fac2,k)*sum2[k]);
    //    std::cout << ccount <<" "<< k<<" "<<sum1[k]<< " "<< sumk1[k]<<" "<<sum2[k] << " "<< sumk2[k]<< std::endl;
    result+=a;
  };

 
   astep = (double) result ;

   //   std::cout <<  astep << std::endl;



   
   };
   
   if (flaglong==0) {ampl2[iev][ib]+= astep;} else Nlong2[iev][ib]++;
   //   std::cout << astep << " ";
     }; // for ib
     //  std::cout << std::endl;
   }; //for iteration
   /* Calculation of the contribution to the amplitude. End. */
  
   if(iev<Nevol)  buffer1->fullevolve(1.0*rapstep);
   //   std::cout << " Next step in rapidity " << std::endl;
    };

  delete buffer1;

  // std::cout <<"**********end event*********" <<std::endl;
  };








    for (int i = 1; i<=Nevol; i++){
      for (int j=0; j<Nb; j++){
	std::cout << ampl2[i][j] << " ";
	ampl[i][j]/=(Nevent-Nlong[i][j]); ampl2[i][j]/=(Nevent-Nlong2[i][j]);
    };
      std::cout << std::endl;
}; 
 



  /*
 for(int i1=0;i1<=Nevol; i1++){
   for(int i2=0; i2<Narray; i2++){
     for(int i3=0; i3<Narray; i3++){
       f1[i1][i2][i3]/=Nevent;
     };
   };
 };
  */


  /* output file */
  std::ofstream of1;
  if1.open(argv[2],std::ios::app);
  if (! if1)
    {
     of1.open(argv[2],std::ios::out);
    };
 of1.precision(15);

 /*
 std::ofstream of1;
 of1.open("eikvsbl.out",std::ios::app);
 if (! of1)
   {
     of1.open("eikvsbl.out",std::ios::out);
   };
 of1.precision(15);

 */
  
 // of1 << "y    ampl_center    Nlong_c    ampl_edge   Nlong_e" << std::endl;

 of1 << "#b, y  ";
 for (int iev=1; iev<=Nevol; iev++){ of1 << rapstart+(iev-1)*rapstep <<" ";}
 of1 << std::endl; 
  for (int j=0; j<Nb; j++){
    of1 << bmin+j*bstep <<" ";
    for (int iev =1; iev<=Nevol; iev++){
    //    of1 <<1+(iev-1)*rapstep <<" "<<ampl[iev] <<" "<<Nlong[iev]<<" "<< ampl2[iev]<<" "<<Nlong2[iev]<<std::endl;
      //    of1 <<1+(iev-1)*rapstep <<" "<< ampl2[iev]<<" "<<Nlong2[iev]<<std::endl;
      of1 << ampl2[iev][j] << " "; 
    }; of1 <<std::endl;
  };
  /*
 for (int in = 0; in<Nevol; in++){
   of1 << in<<" "<<f1[in][Narray/2][Narray/2] << std::endl;
 };
  */
 of1.close();
}


int poisson(double lambda){
  /*lambda = qN*/
  int N = (int) 30*lambda, n=0;
  double q = lambda/N, z;
  for (int i = 0; i<N; i++){z = 1.*rand()/RAND_MAX;
    if (z<q)n++; };
  return n;
};






double makeT(double *T, int n1, int n2)
{
  if ((n1==0)||(n2==0)) return 0;
  double *t; t=T;
  double aT[n1][n2];
      //  std::cout << std::endl;
  for(int i1=0; i1<n1; i1++){
    for(int i2=0; i2<n2; i2++){
      aT[i1][i2]=*t;
      //  std::cout << aT[i1][i2]<<" ";
      t++;
    }; 
      //  std::cout<<std::endl;
};



int count1 = pow(2,n1)-1;
int count2 = pow(2,n2)-1;
//  std::cout << count1<<" "<<count2<<std::endl;                                                                                                           
double prod;
double result=0;
for (unsigned int ig1=1; ig1<=count1; ig1++){
  //  unsigned int ig1 = count2; // code check on computing permanents                                                                                     
  int g1=ig1; // in binary code "g1" is a set of elements s1                                                                                               
  int p[n1],no[n1]; // p --vector of interacting, no -- array with numbers of interacting                                                                  
  int nc1=0; // number of elements in set s1                                                                                                               
  for (int i=0; i<n1; i++){p[i]= g1 &1;
    if (p[i]==1){no[nc1]=i; nc1++;}; g1>>=1;}; // if(nc1<n1-1) no[nc1]=-1;                                                                                 
  double factor = 1;
  double sum[n1]; //nc1 elements of it will be used                                                                                                        
  /* First element of the sum -- s2 = single element #0, doing everything explicitly */
  for(int ii=0; ii<nc1; ii++){ sum[ii]=aT[no[ii]][0]; factor*=sum[ii]; };
  /* s1 contains at least 1 element, so check |s2|<|s1| is not needed */
  factor*=binom(n2-1,nc1-1);
  result+=factor;
  unsigned int g2, g2m = graycode (1);
  unsigned int g3;
  unsigned int gout;
  unsigned int gin;

int nc2=1;
//  for (int iii=0; iii<nc1; iii++)  std::cout << no[iii]<<" ";                                                                                         
//  std::cout << std::endl <<ig1 <<std::endl;                                                                                                           
for (unsigned int ig2=2; ig2<=count2; ig2++ ){
  /*selecting set s2 - gray code*/
  g2 = graycode (ig2);
  g3 = g2 ^ g2m;
  gout = g3 & g2m;
  gin = g3 & g2;
  int flag =0 , nout =-1, nin=-1;
  /* two possibilities -- element added or removed from the sum */
  if(gout>0){//  std::cout << "removed" <<std::endl;                                                                                                      
    do{nout++ ; flag = gout & 1; gout >>=1;}while(flag==0);
    for (int ii=0; ii<nc1; ii++){sum[ii]-=aT[no[ii]][nout];};
    //     std::cout << "removed " <<nout <<std::endl;                                                                                                      
    nc2--;
  }else{//std::cout << "added" <<std::endl;                                                                                                               
    do{nin++ ;
      flag = (gin & 1);
      gin>>=1;
    }while(flag==0);
    //   std::cout << "added "<<nin << std::endl;                                                                                                         
    for (int ii=0; ii<nc1; ii++){ sum[ii]+=aT[no[ii]][nin];};
    nc2++;
  };
  factor=1;
  if(nc2<=nc1){for(int ii=0; ii<nc1; ii++) factor*=sum[ii];
    factor*=(pow(-1,nc2%2+1)*binom(n2-nc2, nc1-nc2));
    result+=factor;}
  g2m = g2;
 };
};
return result;
};





unsigned int graycode(unsigned int g)
{return g^(g>>1);};

unsigned long int binom (unsigned int n, unsigned int k)
{
  unsigned long int factor=1;
  unsigned int an = n;
  for (unsigned int i = 1; i<=k; i++){factor*= an; an--; factor/=i;};
  return factor;
}


/*
double RandGauss()
{
  double r;
  double v1,v2,fac;
        do
          {
            v1 = 2.0 * rand()/RAND_MAX - 1.0;
            v2 = 2.0 * rand()/RAND_MAX - 1.0;
            r = v1*v1 + v2*v2;
          }
        while ( r > 1.0 );

        fac = sqrt(-2.0*log(r)/r);
        return v2*fac;
}
*/
