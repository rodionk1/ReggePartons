#include "pbuffer.hh"
#include "pstack.hh"
#include "parton.hh"

#include "assert.h"


pbuffer::pbuffer()
{  UpperPtr=LowerPtr=0;

  count = 0;
};

pbuffer::~pbuffer()
{
  if (!isEmpty()) //buffer is not empty
    {
      pstack* CurrentPtr=UpperPtr, *TempPtr;
      while(CurrentPtr!=0)
	{ TempPtr=CurrentPtr;
          CurrentPtr=CurrentPtr->NextPtr;
          delete TempPtr;
        };
    };

};

int pbuffer::isEmpty() const
{  return UpperPtr == 0; };

int pbuffer::GetCount()
{return count;};

  void pbuffer::mergestacks(int n1, int n2)
  { // std::cout << "Merging stacks "<<n1<<" and "<< n2 << std::endl; // At fusion the stacks are reordered, so the fusion is always of stack 1 and 2
    if(n2>n1){
      pstack* Temp1; Temp1 = get(n2); 
      pstack* Temp2; Temp2 = get(n1);
      Temp1->LowerPtr->NextPtr = Temp2->UpperPtr;
    // last parton in the first stack is attached to the first parton of the second stack
      Temp2->UpperPtr = 0;
      Temp2->LowerPtr =0;
    // pointers to the partons of the second stack are contained in the first, so we make second buffer empty and delete it
      delete Temp2;
      // now put stack1 with all the partons into the buffer
      put(*Temp1);
      delete Temp1; }else{//another order of buffers
      pstack* Temp1 = get(n1); pstack* Temp2 = get(n2);
      Temp1->LowerPtr->NextPtr = Temp2->UpperPtr;
      // last parton in the first stack is attached to the first parton of the second stack 
      Temp2->UpperPtr = 0;
      Temp2->LowerPtr =0;
      // pointers to the partons of the second stack are contained in the first, so we make second buffer empty and delete it 
      delete Temp2;
      // now put stack1 with all the partons into the buffer    
    put(*Temp1);
      delete Temp1; 
}; 

  };



pstack* pbuffer::read(int n)
{
  assert (n <= count); //check if No does not exceed total number of partons                                                                                 
  // parton with UpperPtr = parton N1, LowerPtr = parton N count                                                                                             
  pstack * TempPtr;
  TempPtr = UpperPtr;
  for (int i=1;i<n; i++) {TempPtr=TempPtr->NextPtr;}; //selecting the right stack                                                                           
  return TempPtr;
};




  pstack * pbuffer::get (int n)
  {
    assert (n <= count); //check if No does not exceed total number of stacks in the buffer
    pstack * TempPtr;
    pstack * TempPtr1;
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
  };


  pstack* pbuffer::getNewstack(const pstack & copy)
  {
    pstack* ptr = new pstack(copy);
    assert(ptr != 0);
    return ptr;
  };



  void pbuffer::put(const pstack & copy)
  {
    count++;
    //    std::cout <<"putting stack into the buffer" <<std::endl;
    pstack* NewPtr= getNewstack(copy);
    //    std::cout <<"created copy of a stack with "<< NewPtr->count<< " partons"<<std::endl;
    if (isEmpty()) // buffer empty
      UpperPtr=LowerPtr=NewPtr;
    else
      { //buffer not empty
	//	std::cout << "buffer not empty, adding stack" <<std::endl;
	NewPtr->NextPtr=UpperPtr;
	UpperPtr=NewPtr;
      };
    //         std::cout << "Now " << count <<" stack(s) in the buffer" <<std::endl;
  };





void pbuffer::diffuse(double t)
{
  if(!isEmpty()){
  pstack* TempPtr = UpperPtr;
  while(TempPtr!=0){
    TempPtr->diffuse(t);
    TempPtr=TempPtr->NextPtr;
    //    std::cout << TempPtr << std::endl;
  };};
};

void pbuffer::evolve(double t)
{
  double tevol = t;
  double tevent,r; 
  while((tevol>0.)&&(!isEmpty())){
    //    std::cout <<"doing evolve with " << count <<" stacks"<<std::endl;
  double lam[count],lamtot=1.e-19;  // lambda and m1 can be different for each stack
  double m1[count],m1tot=1.e-19;
  pstack * TempPtr;
  pstack * TempPtr1;
  double m2, nu;
  if(count>=1){
  TempPtr = UpperPtr;
    m2 = TempPtr->m2;
    nu = TempPtr->nu;
  lam[0]=(TempPtr->lam)*(TempPtr->count); // total splitting probability per unit time, lambda for the current stack
  m1[0]=(TempPtr->m1)*(TempPtr->count); // the same for m1 (the same for all partons in the stack)
  for(int i=2; i<=count; i++)
    {
      TempPtr=TempPtr->NextPtr;
      lam[i-1] = lam[i-2]+(TempPtr->lam)*(TempPtr->count);
      m1[i-1] = m1[i-2]+(TempPtr->m1)*(TempPtr->count);
      //    lamtot+=lam[i-1];
      //  m1tot+=m1[i-1];
    };
  lamtot+=lam[count-1]; // total splittint probability of all partons in all stacks
  m1tot+=m1[count-1];   // total death probability of all partons in all stacks
  };
  //    std::cout <<"lambda and m1 done" <<std::endl;

  /*Creating matrix of interactions within an between stacks*/
  int aint[count][count]; // Works only up to 1770 buffers!
  //  unsigned char aint[count][count];
  //  double aint[count][count];
  double m2nutot=0,m2tot=0,nutot=0;
  for (int i=0;i<count;i++){for (int j=0; j<count; j++){aint[i][j]=0;};};
  TempPtr = UpperPtr; // stack1 , labeled by i
  //  std::cout << TempPtr->count << std::endl;
  for(int i=0; i<count; i++) //perebor po 1 steku
      { TempPtr1 = TempPtr; //stack2, labeled by j
	for (int j=i;j<count;j++){ //perebor po 2 steku
	  //  std::cout <<"making for stacks "<<i+1<<" and "<<j+1<<std::endl;
	  if (i==j) { //esli steki sovpadajut 
	    //  	  std::cout << "steki sovpadajut " <<std::endl;
	    parton* TempPtrPart = TempPtr->UpperPtr; // upper pointer of the stack -- first parton in the stack
	    parton* TempPtrPart1;
	    // std::cout << TempPtr->count << "partons in the stack"<<std::endl;
	    for (int is = 0; is<TempPtr->count;is++){ TempPtrPart1 = TempPtrPart->NextPtr; // another pointer -- next parton of the same stack
	      for (int js = is+1; js<TempPtr->count; js++){		
	  /*computing distance*/
		//		std::cout << is <<" " <<js <<std::endl;
double dist= (TempPtrPart1->x[0]-TempPtrPart->x[0])*(TempPtrPart1->x[0]-TempPtrPart->x[0])+(TempPtrPart1->x[1]-TempPtrPart->x[1])*(TempPtrPart1->x[1]-TempPtrPart->x[1]) ;     
double eps = TempPtr->eps;
 if (dist < eps*eps) { aint[i][j]+=1; m2nutot+=1.;}; 
 TempPtrPart1 = TempPtrPart1->NextPtr;
	      };/*js*/ TempPtrPart=TempPtrPart->NextPtr;}; /*is*/
	    // std::cout << "done" <<std::endl;
	  } else /*if i!=j, interaction between different stacks*/{
parton* TempPtrPart = TempPtr->UpperPtr; // upper pointer of the stack -- first parton in the stack
 parton* TempPtrPart1 ;
 //std::cout << TempPtr->GetCount() <<" " << TempPtr1->GetCount() <<std::endl;
 // std::cout << TempPtrPart->NextPtr << " " << TempPtrPart1->NextPtr << std::endl;
 for (int is = 0; is<TempPtr->count;is++){ TempPtrPart1 = TempPtr1->UpperPtr; // upper pointer of the second stack -- first parton of the second stack
	      for (int js = 0; js<TempPtr1->count; js++){		
	  /*computing distance*/
		//	std::cout << is <<" " <<js <<std::endl;
double dist= (TempPtrPart1->x[0]-TempPtrPart->x[0])*(TempPtrPart1->x[0]-TempPtrPart->x[0])+(TempPtrPart1->x[1]-TempPtrPart->x[1])*(TempPtrPart1->x[1]-TempPtrPart->x[1]) ;    
 double eps = TempPtr->eps; 
if (dist < eps*eps) { aint[i][j]+=1; m2nutot+=1.;}; 
 TempPtrPart1 = TempPtrPart1->NextPtr;
	      };/*is*/ TempPtrPart=TempPtrPart->NextPtr;}; /*js*/  
	  };
        TempPtr1 = TempPtr1->NextPtr;         
	};  
	TempPtr=TempPtr->NextPtr;
      };
    m2tot = m2nutot*m2;
    nutot = m2nutot*nu;
    m2nutot*=(m2+nu);
    /*interaction matrix for stacks completed*/
    /* interaction probabilities per unit time = m2tot + nutot + lamtot +m1tot*/ 

    //    std::cout << "interaction matrices for stacks done " << std::endl;

  /*Defining time before the 1st event*/ 

  r = 1.*rand()/RAND_MAX;
  tevent = -log(r)/((lamtot+m1tot+m2nutot));
  //    std::cout << "lamtot "<<lamtot<<" m1tot "<<m1tot<<" count "<<count<<" r "<<r<<" tevent " <<tevent<<std::endl;
  if (tevent>tevol){
    //  std::cout <<"starting diffusion" <<std::endl;
    diffuse(tevol); tevol = -1.; //no splittng or death   
    //  std::cout << "diffusion" <<std::endl;
  } else   { 
    //    std::cout << "splitting/death/fusion/annihilation" <<std::endl;
      /*if splitting or death - define type of event*/
      r = 1.*rand()/RAND_MAX;
      int type;
      if(r<=lamtot/(lamtot+m1tot+m2nutot)) type = 1; // splitting
      if((r>lamtot/(lamtot+m1tot+m2nutot))&&(r<=(lamtot+m1tot)/(lamtot+m1tot+m2nutot))) type = 2; // death
      if((r>(lamtot+m1tot)/(lamtot+m1tot+m2nutot))&&(r<=(lamtot+m1tot+m2nutot)/(lamtot+m1tot+m2nutot))) type = 3; // fusion or annihilation
 
      int no =0; // number of stack with event / number of parton in the stack
      switch (type){
      case 1:  
	/*diffuse, splitting, define number and splitt*/
	{   diffuse(tevent);
          r=lamtot*rand()/RAND_MAX;
          do{no++;}while(r>lam[no-1]);
	          //              std::cout <<"splitting in stack " <<no <<std::endl;
		  pstack * SplittingStack = get(no); /*extracting stack from buffer*/
		  /*choosing random number to split*/
		  //		  std::cout << "stack extracted " <<std::endl;
		  no = (int) floor(1.*rand()/RAND_MAX*SplittingStack->count)+1;
		  if (no>SplittingStack->count) no = SplittingStack->count; 
		  /* "if" is to prevent a very rare situation with rand()=RAND_MAX which gives no = count+1)*/
		  //		  std::cout << "splitting parton " <<no <<std::endl;
		  SplittingStack->split(no);
		  //		  std::cout << " splitting done, now " << SplittingStack->GetCount() << " partons in the stack "<<std::endl;
		  put(*SplittingStack); /*putting back after splitting*/
		  delete SplittingStack; /*stack is in buffer -- delete to not duplicate*/
		  //  std::cout <<"putback done" <<std::endl;
	  }
	  break;

      case 2:
	/*killing a random parton*/
	{  diffuse(tevent);
        r=m1tot*rand()/RAND_MAX;
        do{no++;}while(r>m1[no-1]);
	//           std::cout << "killing parton from stack " << no << std::endl;                                                                                           
	   /*choosing random parton to kill*/
 	   pstack * KillinStack = get(no); /*extracting stack*/
	   no = (int) floor(1.*rand()/RAND_MAX*KillinStack->count)+1;
	   if ( no > KillinStack->count) no = KillinStack->count;
	   KillinStack->kill(no);
	   if(!(KillinStack->isEmpty())){ put(*KillinStack); /*putting back after killing parton*/};
	   delete KillinStack;
           }

	break;
      case 3:
	{
	  //	  std::cout <<"fusion or annihilation of pair" << std::endl;
	  r=(m2nutot/(nu+m2))*rand()/RAND_MAX;
	  ///	  for(int ii=0; ii<count; ii++){for(int jj=0;jj<count; jj++){std::cout << aint[ii][jj]<<" ";};std::cout<<std::endl;};
	  /* nutot/nu -- number of interactiong pairs */
	  int flagnu =0;
          int inu=-1,jnu=1;
	  long int sumnu =0;
	  //	  double sumnu =0;
          do{inu++; jnu=inu-1; do{ jnu++; sumnu+=aint[inu][jnu]; if(sumnu>r){flagnu=1;};}while((jnu<count-1)&&(flagnu<1));}while(flagnu<1);
	  inu++; jnu++;
	  //	  	  std::cout << "Pair from stacks "<<inu <<" " <<jnu << " out of " << count << "stacks" << std::endl;
	  /*inu, jnu -- numbers of stacks */
	  if (inu==jnu){/*fusion/annihilation of partons from the same stack*/
	    pstack* FusionStack = get(jnu);
	    //    std::cout <<"/*creating interaction matrix for the stack*/" <<jnu<<std::endl;
	    //int aintstack[FusionStack->count][FusionStack->count]; // up to 1770 partons in the stack
	    char aintstack[FusionStack->count][FusionStack->count]; // up to 2500 partons in the stack
	    //    std::cout << FusionStack->count << " partons in the stack" <<std::endl;

	    int npairs = 0; // number of interacting pairs
	    for (int i=0;i<FusionStack->count;i++){for (int j=0; j<FusionStack->count; j++){aintstack[i][j]=0;};};
	
	    parton* TempPart = FusionStack->UpperPtr; parton* TempPart1;
	      //    lam[0]=TempPtr->lam;                                                                                                                               
	      //    m1[0]=TempPtr->m1;                                                                                                                                 
	      for(int i=0; i<FusionStack->count; i++)
		{ TempPart1 = TempPart->NextPtr;
		  for (int j=i+1;j<FusionStack->count;j++){
		    /*computing distance*/
         double dist= (TempPart1->x[0]-TempPart->x[0])*(TempPart1->x[0]-TempPart->x[0])+(TempPart1->x[1]-TempPart->x[1])*(TempPart1->x[1]-TempPart->x[1]);
	 if (dist < FusionStack->eps*FusionStack->eps) { aintstack[i][j]=1; npairs++;} else{aintstack[i][j]=0;};
		    TempPart1 = TempPart1->NextPtr;
		  };
		  TempPart=TempPart->NextPtr;
		};
	      //   std::cout <<"/*defined number of interactions and interaction matrix*/" <<std::endl;
	      diffuse(t); // diffuse the buffer and the stack
	      //   std::cout <<"diffusion of buffer done" <<std::endl;
	      FusionStack->diffuse(t);
	      //   std::cout <<"diffusion of stack done" <<std::endl;
	      /*now define pair of partons to annihilate or to fuse*/ 
	      
	      r=npairs*(1.*rand())/RAND_MAX;
	      int flagnu =0;
	      int inus=-1,jnus=1;
	      int sumnus =0;
	      //	      double sumnus =0;
	      //    std::cout <<"Selecting pair" << std::endl;
	      do{inus++; jnus=inus; do{jnus++; sumnus+=aintstack[inus][jnus];// std::cout << inus <<" " <<jnus << std::endl; 
                if(sumnus>r){flagnu=1;};}while((jnus<FusionStack->count-1)&&(flagnu<1));}while(flagnu<1);
	      //max value of jnu = count-1, number of partons is inu/ju +1
	      inus++; jnus++;
	      //  std::cout <<"fusion of pair " << inus <<" "<< jnus <<" of "<<FusionStack->count <<" partons"<<std::endl;                                     
	      //     for (int ii=0;ii<count;ii++){for(int jj=0;jj<count;jj++){std::cout<<aint[ii][jj]<<" ";};std::cout<<std::endl;};                                  
	      parton* parton1 = FusionStack->get(inus);
	      /*i<j always, so j must be reduced by 1 after we took away parton #i*/
	      jnus--;
	      parton* parton2=FusionStack->get(jnus);
	      r = (FusionStack->m2+FusionStack->nu)*rand()/RAND_MAX;
	      if(r<FusionStack->nu){//fusion
		//		std::cout << "fusion"<<std::endl;
	      double x=0.5*(parton1->getX(1)+parton2->getX(1)), y=0.5*(parton1->getX(2)+parton2->getX(2));
	      /*arranging parton after fusion in between the two which we remove*/
	      //      std::cout<<parton1->getX(1)<<" "<<parton1->getX(2) <<std::endl;                                                                          
	      //      std::cout<<parton2->getX(1)<<" "<<parton2->getX(2) <<std::endl;                                                                          
	      delete parton1;
	      parton1 = new parton(x,y);
	      delete parton2;
	      FusionStack->put(*parton1);
	      delete parton1;
	      put(*FusionStack);/*putting fusion stack into the buffer*/
	      delete FusionStack;}else{ //annihilation of a pair
		//		std::cout << "same stack annihilation" <<std::endl;
		delete parton1; delete parton2;
		if(!FusionStack->isEmpty()) put(*FusionStack);
		delete FusionStack;
	      };
	  }else{  /*inu!=jnu*/ /* fusion/annihilation of partons from different stacks*/
	    //  std::cout <<"fusion from different stacks" <<std::endl;
	    //	    std::cout << inu <<" " << jnu << " " <<count <<std::endl;
	    pstack* FStack1=get(jnu); // jnu > inu so taking it first - otherwise need get(inu) jnu-- get(jnu)
	    pstack* FStack2=get(inu);
	    int count1 = FStack1->count;
	    int count2 = FStack2->count;
	    char aintstack[count1][count2]; // char -- to allow up to 2500 partons
	    int nint=0;
	    for (int is = 0; is <count1; is++){for(int js = 0; js<count2 ; js++){
		aintstack[is][js]=0;};};
	    parton* PartPtr1 = FStack1->UpperPtr;
	    parton* PartPtr2;
	    /*creating matrix of interstack interactions*/
	    
	    for(int is = 0; is <count1; is++){
	    PartPtr2 = FStack2->UpperPtr;
	      for(int js = 0; js<count2 ; js++){
		double dist = pow(PartPtr1->x[0]-PartPtr2->x[0],2)+pow(PartPtr1->x[1]-PartPtr2->x[1],2);
		double eps = FStack1->eps;
		if (dist<eps*eps) {aintstack[is][js]=1; nint++;} else aintstack[is][js]=0;
		PartPtr2=PartPtr2->NextPtr;}; 
	    PartPtr1=PartPtr1->NextPtr;};
	    //	    std::cout << "matrix of interactions done" <<std::endl;
	    r = (1.*rand()/RAND_MAX)*nint;
	    int inus = -1, jnus;
	    double sumnus = 0;
	    do{inus++; jnus=-1; do{jnus++; sumnus+=aintstack[inus][jnus]; if(sumnus>r){flagnu=1;};}while((jnus<FStack2->count-1)&&(flagnu<1));}while(flagnu<1);
	      //max value of jnu = count-1, number of partons is inu/ju +1
	      inus++; jnus++;
	      /*fusing / annihilating partons inu, jnu from stacks 1 and 2*/
	      
	      diffuse(t);
	      FStack1->diffuse(t);
	      FStack2->diffuse(t);
	      r=(1.*rand()/RAND_MAX)*(m2+nu);
	      if(r<m2){/*annihilation*/
		//	std::cout <<"annihilation from different stacks" <<std::endl;	
      		FStack1->kill(inus);
		FStack2->kill(jnus);
		if(!(FStack1->isEmpty())) put(*FStack1);
		if(!(FStack2->isEmpty())) put(*FStack2);
		delete FStack1; delete FStack2;
	      }else{/*fusion*/
		parton* parton1 = FStack1->get(inus);
		parton* parton2 = FStack2->get(jnus);
		double x = 0.5*(parton1->x[0]+parton2->x[0]);
		double y = 0.5*(parton1->x[1]+parton2->x[1]);
		delete parton1;
		delete parton2;
		parton1 = new parton(x,y);
		FStack1->put(*parton1);
		put(*FStack1);
		if(!(FStack2->isEmpty())){ put(*FStack2); mergestacks(1,2);};  /*buffer is FIFO, so the stacks are #1 and #2*/
		delete FStack1;
		delete FStack2;
	      };
	      
	  };
	}
      break;
      };
  tevol-=tevent;
  };
  //    std::cout <<"left for evolution "<<tevol <<" t "<<t <<std::endl;
  };
  //  std::cout <<" end of evolution" <<std::endl;
}




void pbuffer::fullevolve(double t){
  double tevol=t;
  while (tevol>0){
    //    std::cout <<"*"<<std::endl;                                                                                                                                              
    if (count>=1)
      {double tmax = (0.2*UpperPtr->eps*UpperPtr->eps)/(UpperPtr->alp);
	//      std::cout<<tmax<<std::endl;                                                                                                                                          
	if(tmax>tevol) {
	  evolve(tevol);
	  //      std::cout << "evolving by "<<tevol<<std::endl;                                                                                                                     
	  tevol=-1;}else{
	  // std::cout <<"*"<<std::endl;                                                                                                                                             
	  evolve(tmax);
	  //    std::cout <<"evolving by "<<tmax<<std::endl;                                                                                                                             
	  tevol-=tmax;
	  // std::cout << tevol<<" left for evolution" << std::endl;
	}
      }else{tevol=-1;};
  }; //while tevol>0;               
};

