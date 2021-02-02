#include "stdlib.h"
#include "time.h"
#include "math.h"



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
