#include "parton.hh"
#include "assert.h"
#include <iostream>
//#include "RGauss.hh"

//DESTRUCTOR                                                                                                                                                 
//double RandGauss();

parton::~parton()
{
};


//CONSTRUCTOR, RELOADED                                                                                                                                      

//zero constructor     
/*
parton::parton(double x1, double x2, double *param)
{NextPtr=0;
  x[0] = x1;
  x[1] = x2;
  double * pparam = param;
  alp = *pparam; 
  alp*=2.; //diffusion coefficient for partons is 2x larger than alpha' to reproduce regge behaviour with alpha'
  pparam++;
  lam = *pparam;
  pparam++;
  m1 = *pparam;
  //  beta = 0;
  //  ltime = 0;
  //  tmin = 0.;
}
*/

parton::parton(double x1, double x2)
{NextPtr=0;
  x[0] = x1;
  x[1] = x2;
  //  beta = 0;
  //  ltime = 0;
  //  tmin = 0.;
};



//copy constructor                                                                                                                                           
parton::parton(const parton & copy)
{

  //  NextPtr=copy.NextPtr;
  NextPtr=0;
  x[0] = copy.x[0];
  x[1] = copy.x[1];
  /*  alp = copy.alp;
  lam = copy.lam;
  m1 = copy.m1;*/

};

/*
void parton::diffuse(double t)
{
  x[0]+=RandGauss()*sqrt(alp*t);
  x[1]+=RandGauss()*sqrt(alp*t);
};
*/

double parton::getX(int i)
{
  assert(i<3);
  return x[i-1];
};

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
};
*/
