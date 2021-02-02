

#ifndef PSTACK_H
#define PSTACK_H

//#include "pbuffer.hh"

class pbuffer;
class parton;

class pstack
{
  friend class pbuffer;
  //  static const double m2=0.0; // parton pair death probability
  //  static const double nu=0.0; // parton fusion probability
  //  static const double eps=0.02; // parton interaction distance

public:
  pstack(double *params); //constructor
  pstack(const pstack & copy); //copy constructor
  ~pstack(); // destructor                                                                                                                                 
  void put(const parton &);
  parton* get();
  parton* get(int n);
  parton* read(int n); //read parton number n from the stack
  parton* getNewparton(const parton &);
  void diffuse(double t); // diffusion of whole stack of partons
  void split(int); // split parton
  void kill(int); // kill parton
  //  void evolve (double t); // evolve partonic distributions
  //  void fullevolve(double t); // evolve distributions by small times paying attention to the fusion of partons
  int isEmpty() const;
  int GetCount(); //get number of partons in the stack
  double GetEps();
private:
  double m2,nu,eps,alp,lam,m1;
  int count; // number of partons in the stack
  parton* UpperPtr; // pointer to the last added parton in the stack
  parton* LowerPtr; // pointer to the first added parton                                                                        
  //  parton* getNewparton(const parton &);
  pstack* NextPtr; // pointer to the new stack in the buffer
};


#endif
