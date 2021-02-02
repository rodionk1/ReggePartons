
#ifndef PBUFFER_H
#define PBUFFER_H


class pstack;


class pbuffer
{
public:
  pbuffer();
  ~pbuffer();
  void mergestacks(int n1, int n2);
  pstack* get(int n);
  void put(const pstack &);
  pstack* read(int n); //read stack number n from the buffer
  pstack* getNewstack(const pstack &);
  void diffuse(double t); // diffusion of whole buffer of partons                   
  void evolve (double t); // evolve partonic distributions
  void fullevolve(double t); // evolve distributions by small times paying attention to the fusion of partons                
  int GetCount(); // number of stacks in the buffer
  int isEmpty() const;
private:
  int count; //number of stacks in the buffer
  pstack* UpperPtr;
  pstack* LowerPtr;
};

#endif
