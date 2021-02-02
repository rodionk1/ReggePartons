
#include <fstream>
#include <iostream>
#ifndef PARTON_H
#define PARTON_H
//#include "time.h"
//#include "math.h"
//#include "pstack.hh"
//#include "pbuffer.hh"
class pstack;
class pbuffer;

class parton
{
  friend class pstack;
  friend class pbuffer;
  //  static const double alp=0.0485; // diffusion coefficient \alpha'
  //  static const double lam=0.0; // splitting probability
  //  static const double m1=0.0;  // death probability
  //  static const double m2=0.1;  // pairwise death probability
  //  static const double nu=0.1;  // fusion probability
  //  static const double eps=0.1; // parton interaction radius 
  //  double alp,lam,m1; //mmoved all the parameters to the definition of a stack
public :
  //  parton(double x, double y, double *param);//construction with coordinates
  parton(double x, double y);
  parton(const parton & copy); // copy constructor  
 ~parton(); // destructor
  // void diffuse(double t); // diffusion in b-plane for time t  --> moved to evolution of a stack
  double getX(int i); // get coordinate number i at given time
private : 
  parton * NextPtr;
  double x[2]; // coordinates
};

#endif
