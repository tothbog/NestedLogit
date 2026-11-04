#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "GORoutines.h"

extern int    n;          // number of demand points
extern REAL   b1, b2;     // parameters of the distance function
extern MATRIX p;          // locations of the demand points
extern VECTOR w;          // buying power at the demand points
extern REAL   delta;      // the forbidden radius around demand point i is w(i)/delta

/*----------------------------------------------------------------------------*/
void InitBox ( PBOX & pb, const INTERVAL_VECTOR & XLim, HTscalar_FctPtr ff,
               int & NEvalF )
{
 pb->X = XLim;
 fEvalH(ff, pb->X, pb->FX); NEvalF++;
}

/*----------------------------------------------------------------------------*/
void InitBoxGrad ( PBOX & pb, INTERVAL_VECTOR & GX, const INTERVAL_VECTOR & XLim,
		HTscalar_FctPtr ff, int & NEvalG )
{
 pb->X = XLim;
 fgEvalH(ff, pb->X, pb->FX, GX); NEvalG++;
}

/*----------------------------------------------------------------------------*/
void BaumannPoint(const PBOX & pb, INTERVAL_VECTOR & GX, PBOX & pc)
{
 int Dim = Dimension(*pb);
 VECTOR wGX(Dim), cBaum(Dim);

 wGX = Diam(GX);

 for(int j=1;j<=Dim;j++)
    {
     if (0.0 < Inf(GX(j)))  cBaum(j) = Inf(pb->X(j));
     else
       if (0.0 > Sup(GX(j))) cBaum(j) = Sup(pb->X(j));
       else
         cBaum(j) = ( Sup(GX(j))*Inf(pb->X(j)) - Inf(GX(j))*Sup(pb->X(j)) ) / wGX(j);
    }
 pc->X = Hull(cBaum);
}

/*----------------------------------------------------------------------------*/
bool UpdateBetterUbf (const PBOX & pb, REAL & BetterUbf, VECTOR & BetterPoint, 
                      HTscalar_FctPtr ff, int & NEvalFx, PBOX & pc)
{
 if ( BetterUbf > Sup(pc->FX) )
    {
     BetterUbf=Sup(pc->FX);
     BetterPoint=Mid(pc->X); 

     cerr << "BetterUbf = " << BetterUbf << endl;
     cerr << "\tat BetterPoint = " << BetterPoint << endl;
     return true;
    } 

 return false;
}

/*----------------------------------------------------------------------------*/
void CutOffTest (PBLTREE pblt, REAL BetterUbf, int & NDelCTT,
                 const INTERVAL_VECTOR & XLim, INT WWidth, 
                 bool Draw, PCHAR color )
{
 int NDel;

 NDel=pblt->Dimension();
 
 if ( BetterUbf < pblt->GetKeyCut() )
    pblt->DelGTKey (BetterUbf, XLim, WWidth, Draw, color);

 NDel-=pblt->Dimension();
 NDelCTT+=NDel;	    
}

/*----------------------------------------------------------------------------*/
bool CutOffTheActualBox ( PBOX & pb, REAL BetterUbf, int & NDelMPT,
                    const INTERVAL_VECTOR & XLim, INT WWidth, bool Draw)
{
 bool COTest=false;
 
 if ( BetterUbf < Inf(pb->FX) )
    {
     if (Draw)
        pb->Draw(XLim,WWidth,"Gray60");
     delete pb;
     pb=NULL;
     NDelMPT++;
     COTest=true;
    } 
    
 return COTest;
}

/*----------------------------------------------------------------------------*/
bool MonoTestFeas(PBOX & pb, INTERVAL_VECTOR & GX, int & NDelMT, const INTERVAL_VECTOR & XLim, 
			INT WWidth, bool Draw)
{
 bool Monotone=false;
 int NDim = Dimension(pb->X);
 INTERVAL_VECTOR Y(NDim);  // the box restricted to the border facet
 PBOX pt=NULL;
 pt = new BOX(*pb);        // the facet where the minimum over the box is attained
 Y=pb->X;
 bool border = false;
 
 pb->M=false;
 for (int i=1; i<=NDim; i++ )
     if (  !((REAL)0.0 <= GX(i))  )   // the gradient does not contain 0
	{
	   pb->M = true;
	   if ( Inf(GX(i)) > 0.0 )         // increasing: the minimum is at the lower end
		{
		 pt->X(i) = Inf(pb->X(i));
		 if (Inf(pb->X(i)) <= Inf(XLim(i)))  // the box is on the border
		    {
		     Y(i)= Hull(Inf(XLim(i)),Inf(pb->X(i)));
		     border = true;
		    }
		}
	   else                            // decreasing: the minimum is at the upper end
		{
		 pt->X(i) = Sup(pb->X(i));
		 if (Sup(pb->X(i)) >= Sup(XLim(i)) )  // the box is on the border
		   {
		    Y(i) = Hull(Sup(XLim(i)),Sup(pb->X(i)));
		    border = true;
		   }
		}
	} 

 // monotone and the facet with the minimum is feasible
 if ( pb->M && ( pb->feas || Feasible(pt)==1) ) 
    {
     if (!border)
       { 
	if (Draw)
	  {
	   pb->UnDraw();
	   pb->Draw(XLim,WWidth,"Orange");
	  } 
	delete pb;
	pb=NULL;
	Monotone=true;
	NDelMT++;
       }
     else pb->X = Y;
    } 
 
 delete pt; pt = NULL;
 return Monotone;
}

/*----------------------------------------------------------------------------*/
/* The new facility must be outside the circle of radius w(i)/delta around     */
/* every demand point i: b1*(x1-p_i1)^2 + b2*(x2-p_i2)^2 >= (w(i)/delta)^2.    */
/* fe(i) records the constraints already known to hold on the whole box.       */
/*----------------------------------------------------------------------------*/
int Feasible (PBOX & pb)
{
  if (pb->feas==true)
  return 1;

  INTERVAL d, dd, tmp;
  bool feasible=true;
  for (int i=1; i <= n; i++)
    if ( !pb->fe(i) )
    {
      tmp = pb->X(1) - p(i,1);
      d = Hull(pow(Inf(tmp), 2), pow(Sup(tmp), 2));
      if ( 0 <= tmp )
         d = Hull(0.0, d);
      dd = b1 * d;
      tmp = pb->X(2) - p(i,2);
      d = Hull(pow(Inf(tmp), 2), pow(Sup(tmp), 2));
      if ( 0 <= tmp )
         d = Hull(0.0, d);	
      
      dd = dd + b2 * d - pow(w(i)/delta, 2);
      
      if ( Sup(dd) < 0.0 )        // the box is infeasible, can be deleted
         return -1;
      else if ( Inf(dd) >= 0.0 ) 
        pb->fe(i) = true;
      else
        feasible = false;
    }

  if (feasible)
  {
    pb->feas=true;
    return 1;
  }
  else
    return 0;
}
