#include "IntervalVector.h"
#include "IntervalVectorPlusz.h"

/*----------------------------------------------------------------------------*/

REAL   SumDiam   (const INTERVAL_VECTOR & X)
{
  REAL sum=0.0;
  for (int i=1; i<=Dimension(X); i++)
  	sum += Diam(X(i));
	
  return sum;
}

/*----------------------------------------------------------------------------*/
REAL   MaxDiam   (const INTERVAL_VECTOR & X)
{
  REAL max=0.0;
  for (int i=1; i<=Dimension(X); i++)
  	if ( Diam(X(i)) > max )
	    max = Diam(X(i));
	
  return max;
}

/*----------------------------------------------------------------------------*/

INT    NDiamGtEp (const INTERVAL_VECTOR & a, REAL Eps)
{
  if (MaxDiam(a) > Eps) return 1; else return 0;
}

