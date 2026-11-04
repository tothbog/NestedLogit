#include <stdlib.h>
#include "Sorting.h"

/*--------------------------------------------------------------------------*/
/*Order porder from upper to lower values of pdist.                         */
/*porder must be initialized to natural order.                              */
/*--------------------------------------------------------------------------*/
void Bubble(PINT porder, PREAL pdist, INT NDim)
{
 INT i,j,Temp;
 
 for (i=0;i<NDim-1;i++)
     for (j=i+1;j<NDim;j++)
         {
	  if (pdist[porder[i]] < pdist[porder[j]])
	     {
	      Temp=porder[i];
              porder[i]=porder[j];
              porder[j]=Temp;
	     }
	 }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
