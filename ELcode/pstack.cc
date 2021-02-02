#include "parton.hh"
#include "pstack.hh"
#include "assert.h"
//#include "pbuffer.hh"



double pstack::GetEps()
{return eps;
};

int pstack::GetCount()
{ return count;
};


/* constructor */
pstack::pstack( double * params)
{ double * param = params;
  m2 = *param; param++;
  nu = *param; param++;
  eps = *param;param++;
  alp = *param; alp*=2.;   param++;
  lam = *param; param++;
  m1 = *param; param++;
  UpperPtr=LowerPtr=0;
  NextPtr=0;
  count = 0;
};


pstack::pstack(const pstack& copy)
{
  m2=copy.m2;
  nu=copy.nu;
  eps=copy.eps;
  alp=copy.alp;
  lam=copy.lam;
  m1=copy.m1;
  UpperPtr=LowerPtr=0;
  NextPtr=0;
  count =0;
  //  std::cout << copy.count <<std::endl;
  parton* NewPart = copy.UpperPtr;
  while (NewPart!=0){
  //  std::cout << "creating new stack " <<count << std::endl;
    put(*NewPart);
  NewPart = NewPart->NextPtr;
  };
};


/* destructor */
pstack::~pstack()
{
  if (!isEmpty()) //stack not empty                                                                                                                         
    {
      parton* CurrentPtr=UpperPtr, *TempPtr;
      while(CurrentPtr!=0)
	{ TempPtr=CurrentPtr;
	  CurrentPtr=CurrentPtr->NextPtr;
	  delete TempPtr;
	};
    };
  //  std::cout << "No more strings" << std::endl;
};


void pstack::put(const parton & copy)
{
  count++;
  parton* NewPtr= getNewparton(copy);
  if (isEmpty()) // stack empty                                                                                                                             
    UpperPtr=LowerPtr=NewPtr;
  else
    { //stack not empty                                                                                                                                     
      NewPtr->NextPtr=UpperPtr;
      UpperPtr=NewPtr;
    };
  //  std::cout << "Now " << count <<" parton(s) in the stack" <<std::endl;
};


parton* pstack::get()
{
  if (isEmpty()) return 0;
  else
    {
      parton* TempPtr=UpperPtr;
      UpperPtr=TempPtr->NextPtr; // reducing the buffer                                                                                                         
      count--;
      //      std::cout << count << std::endl;                                                                                                                     
      return TempPtr; // and returning the element which we subtracted                                                                                          
    };
};


int pstack::isEmpty() const
{
  return UpperPtr == 0;
}


parton* pstack::getNewparton(const parton & copy)
{
  parton* ptr = new parton(copy);
  assert(ptr != 0);
  return ptr;
};

void pstack::diffuse(double t)
{
  if(!isEmpty()){
  parton * TempPtr; 
  TempPtr = UpperPtr;
  // Diffusing partons one by one:
  //for (int i=1;i<=count; i++)
  while (TempPtr!=0){TempPtr->x[0]+=RandGauss()*sqrt(alp*t);
    TempPtr->x[1]+=RandGauss()*sqrt(alp*t);
    TempPtr=TempPtr->NextPtr;
    //    std::cout<<count <<" partons "<<TempPtr<<std::endl;
};};
};

void pstack::split(int n)
//split parton No n
{
  assert (n <= count); //check if No does not exceed total number of partons
  // parton with UpperPtr = parton N1, LowerPtr = parton N count
  parton * TempPtr;
  TempPtr = UpperPtr;
  for (int i=1;i<n; i++) {TempPtr=TempPtr->NextPtr;}; //selecting the right parton
  put(*TempPtr);//adding another copy of parton to the stack = 'splitting'
  /*count is increased in put()*/ 
}

parton* pstack::read(int n)
{
  assert (n <= count); //check if No does not exceed total number of partons                                                                 
  // parton with UpperPtr = parton N1, LowerPtr = parton N count                                                                             
  parton * TempPtr;
  TempPtr = UpperPtr;
  for (int i=1;i<n; i++) {TempPtr=TempPtr->NextPtr;}; //selecting the right parton
  return TempPtr;
}


parton * pstack::get (int n)
{
  assert (n <= count); //check if No does not exceed total number of partons 
  parton * TempPtr;
  parton * TempPtr1;
  //  std::cout << n <<" " << count <<std::endl;
  //  std::cout << UpperPtr << " " << LowerPtr << std::endl;                               
  TempPtr = UpperPtr;
  for (int i=1;i<n; i++) {TempPtr1=TempPtr; TempPtr=TempPtr->NextPtr;
    // std::cout << i <<" "<< TempPtr1 <<" " <<TempPtr1->NextPtr<<" "<<TempPtr << std::endl; 
  }; //selecting the right parton                                                                                       
  //  std::cout << TempPtr1 <<" " << TempPtr << std::endl;                                                             
  if(n==1){
    UpperPtr = TempPtr->NextPtr;}
  else if (n==count){
    LowerPtr=TempPtr1;
    //   std::cout << "-->"<<LowerPtr << std::endl;                                                               
    TempPtr1->NextPtr=0;
  } else {
    TempPtr1->NextPtr = TempPtr->NextPtr;}
  count--;
  return TempPtr;
  //  std::cout << UpperPtr << " " << LowerPtr << std::endl;                                                                                                
}

void pstack::kill(int n)
//kill parton No n
{
  assert (n <= count); //check if No does not exceed total number of partons
  parton * TempPtr;
  parton * TempPtr1;
  //  std::cout << n <<" " << count <<std::endl;
  //  std::cout << UpperPtr << " " << LowerPtr << std::endl;
  TempPtr = UpperPtr;
  for (int i=1;i<n; i++) {TempPtr1=TempPtr; TempPtr=TempPtr->NextPtr;
    // std::cout << i <<" "<< TempPtr1 <<" " <<TempPtr1->NextPtr<<" "<<TempPtr << std::endl;
}; //selecting the right parton
  //  std::cout << TempPtr1 <<" " << TempPtr << std::endl; 
  if(n==1){
    UpperPtr = TempPtr->NextPtr; delete TempPtr;}
  else if (n==count){
   LowerPtr=TempPtr1; 
   //   std::cout << "-->"<<LowerPtr << std::endl; 
 TempPtr1->NextPtr=0; delete TempPtr;
  } else {
    TempPtr1->NextPtr = TempPtr->NextPtr; delete TempPtr;}
  count--;
  //  std::cout << UpperPtr << " " << LowerPtr << std::endl;
}

/*************
void pstack::evolve(double t)
{
  double tevol = t;
  double tevent,r; 
 while(tevol>0.){
  double lam[count],lamtot=1.e-19;
  double m1[count],m1tot=1.e-19;
  parton * TempPtr;
  parton * TempPtr1;
  if(count>=1){
  TempPtr = UpperPtr;
  for(int i=2; i<=count; i++)
    {
      lam[i-1] = lam[i-2]+TempPtr->lam;
      m1[i-1] = m1[i-2]+TempPtr->m1;
      TempPtr=TempPtr->NextPtr;
      //    lamtot+=lam[i-1];
      //  m1tot+=m1[i-1];
    };
  lamtot=lam[count-1];
  m1tot=m1[count-1];
  };
  //  std::cout <<"lambda and m1 done" <<std::endl;
  **********************/
  /*Creating matrix of interactions*/
  /**************************
  double aint[count][count];
  double m2nutot=0,m2tot=0,nutot=0;
  for (int i=0;i<count;i++){for (int j=0; j<count; j++){aint[i][j]=0;};};
  if(count>=2){
    TempPtr = UpperPtr;
    //    lam[0]=TempPtr->lam;
    //    m1[0]=TempPtr->m1;
    for(int i=0; i<count; i++)
      { TempPtr1 = TempPtr->NextPtr;
	for (int j=i+1;j<count;j++){
	******************/
	  /*computing distance*/
  /***************8
double dist= (TempPtr1->getX(1)-TempPtr->getX(1))*(TempPtr1->getX(1)-TempPtr->getX(1))+(TempPtr1->getX(2)-TempPtr->getX(2))*(TempPtr1->getX(2)-TempPtr->getX(2)) ;    
 if (dist < eps) { aint[i][j]=1; m2nutot+=1.;}; 
       TempPtr1 = TempPtr1->NextPtr;         
	};  
	TempPtr=TempPtr->NextPtr;
      };
    m2tot = m2nutot*m2;
    nutot = m2nutot*nu;
    m2nutot*=(m2+nu);
  };
  *****************/
  /*Defining time before the 1st event*/ 
  /***********************8
  r = 1.*rand()/RAND_MAX;
  tevent = -log(r)/((lamtot+m1tot+m2nutot));
  //   std::cout << "lamtot "<<lamtot<<" m1tot "<<m1tot<<" count "<<count<<" r "<<r<<" tevent " <<tevent<<std::endl;
  if (tevent>tevol){
    diffuse(tevol); tevol = -1.; //no splittng or death   
  } else   { 
      diffuse(tevent);
  *****************/
      /*if splitting or death - define type of event*/
  /*********************8
r = 1.*rand()/RAND_MAX;
      int type;
      if(r<=lamtot/(lamtot+m1tot+m2nutot)) type = 1; // splitting
      if((r>lamtot/(lamtot+m1tot+m2nutot))&&(r<=(lamtot+m1tot)/(lamtot+m1tot+m2nutot))) type = 2; // death
      if((r>(lamtot+m1tot)/(lamtot+m1tot+m2nutot))&&(r<=(lamtot+m1tot+nutot)/(lamtot+m1tot+m2nutot))) type = 3; // fusion
      if(r>(lamtot+m1tot+nutot)/(lamtot+m1tot+m2nutot)) type = 4; // annihilation
 
      int no =0;
      switch (type){
      case 1:  *******************/
	/*splitting, define number and splitt*/
  /*********************************8
	{
          r=lamtot*rand()/RAND_MAX;
          do{no++;}while(r>lam[no-1]);
	  //	  std::cout <<"splitting parton " <<no <<std::endl;
	  split(no);}
	  break;

      case 2:
  ****************/
	/*killing a random parton*/
  /******************************
	{
        r=m1tot*rand()/RAND_MAX;
        do{no++;}while(r>m1[no-1]);
        //      std::cout << "killing parton " << no << std::endl;                                                                                           
        kill(no);}

	break;
      case 3:
	{
	  //  std::cout <<"fusion of pair" << std::endl;
	  r=(nutot/nu)*rand()/RAND_MAX;
	  **************************/
	  /* nutot/nu -- number of interactiong pairs */
  /***********************************
	  int flagnu =0;
          int inu=-1,jnu=1;
	  double sumnu =0;
          do{inu++; jnu=inu; do{jnu++; sumnu+=aint[inu][jnu]; if(sumnu>r){flagnu=1;};}while((jnu<count)&&(flagnu<1));}while(flagnu<1);
	  inu++; jnu++;
     //	  std::cout <<"fusion of pair " << inu <<" "<< jnu <<" of "<<count <<" partons"<<std::endl;
     //	    for (int ii=0;ii<count;ii++){for(int jj=0;jj<count;jj++){std::cout<<aint[ii][jj]<<" ";};std::cout<<std::endl;};
	    parton* parton1 = get(inu);
	    ******************/
	    /*i<j always, so j must be reduced by 1 after we took away parton #i*/
  /***********************
 jnu--;
	    parton* parton2=get(jnu);
            double x=0.5*(parton1->getX(1)+parton2->getX(1)), y=0.5*(parton1->getX(2)+parton2->getX(2));
  ******************/
    /*arranging parton after fusion in between the two which we remove*/
	    //	    std::cout<<parton1->getX(1)<<" "<<parton1->getX(2) <<std::endl;
	    //	    std::cout<<parton2->getX(1)<<" "<<parton2->getX(2) <<std::endl;
  /**************************************	 
   delete parton1;
	    double newpar[3];
            newpar[0]=parton2->alp; newpar[1]=parton2->lam; newpar[2]=parton2->m1;
	    parton1 = new parton(x,y, &newpar[0]);
	    delete parton2;
	    put(*parton1);
	}
	    break;

      case 4:
  ****************************/
	/* annihilation */
/***********8
	{	r=(m2tot/m2)*rand()/RAND_MAX;
******************/
	/*m2tot/m2 - total number of interacting pairs*/
  /************************
	int flagm2 =0;
	int im2=-1,jm2=1;
	double summ2 =0;
	do{im2++; jm2=im2; do{jm2++; summ2+=aint[im2][jm2]; if(summ2>r){flagm2=1;};}while((jm2<count)&&(flagm2<1));}while(flagm2<1);
	im2++; jm2++;
	//	std::cout <<"annihilation of pair " << im2 <<" "<< jm2 <<" of "<<count <<" partons"<<std::endl;
	    //	    std::cout << " killing " << i <<std::endl;
            kill(im2);	    ****************************/  /*i<j always, so j must be reduced*/ 
  /************************
 jm2--; 
	    //	    std::cout << "killing "<< j <<std::endl; 
	    kill(jm2);}
      break;
      };
  tevol-=tevent;
  };
  //  std::cout <<"left for evolution "<<tevol <<" t "<<t <<std::endl;
 };
 // std::cout <<" end of evolution" <<std::endl;
}
  *************************/
  /**********************
void pstack::fullevolve (double t)
{
  double tevol=t;
  while (tevol>0){
    //    std::cout <<"*"<<std::endl;
  if (count>=1) 
    {double tmax = 0.05/sqrt(UpperPtr->alp);
      //      std::cout<<tmax<<std::endl;
      if(tmax>t) {tevol = t;    
     	evolve(tevol);
	//	std::cout << "evolving by "<<tevol<<std::endl; 
	tevol=-1;}else{
	// std::cout <<"*"<<std::endl;
    evolve(tmax); 
    //    std::cout <<"evolving by "<<tmax<<std::endl; 
    tevol-=tmax;}
    }else{tevol=-1;};
  }; //while tevol>0;
}
*************************/

