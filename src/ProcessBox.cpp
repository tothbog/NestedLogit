#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "Error.h"
#include "Sorting.h"
#include "ProcessBox.h"

/*----------------------------------------------------------------------------*/
void CheckBox (const PBLTREE & pbltwork, const PBLTREE & pbltfinal, 
		PBOX & pb, INTERVAL_VECTOR & GX, PBOX & pcF,
                HTscalar_FctPtr f, REAL & BetterUbf, VECTOR & BetterPoint,
		REAL Epsilon,
                const INTERVAL_VECTOR & XLim, INT WWidth, bool Draw,
		bool MTest, bool EFatherGX, int & NEvalFx, int & NEvalFX, 
		int & NEvalG, 
		int & NDelMPT, int & NDelMT, int & NDelCTT, int & NDelFathGX)
{
 REAL Pf;
 PBOX pc=NULL;  // the box of the Baumann point
 INTERVAL FX;

 if (MTest && EFatherGX)
    {
     // Bound the subbox with the centered form of the father box (no new evaluation).
     Intersection(pb->FX, pb->FX, (pcF->FX + (GX)*(pb->X - pcF->X)));    
     if (CutOffTheActualBox ( pb, BetterUbf, NDelMPT, XLim, WWidth, Draw))
        {
         NDelFathGX++;
         return; 
        }
    }
 else
    {
     // Interval evaluation of the objective; the father's inclusion may be tighter.
     fEvalH(f, pb->X, FX); NEvalFX++;
     Intersection(pb->FX, pb->FX, FX);
    }

 // Cut-off test for the actual box.
 if ( CutOffTheActualBox ( pb, BetterUbf, NDelMPT, XLim, WWidth, Draw)) 
    return;

 // Feasibility test.
 if ( Feasible(pb) == -1 )
    {
     delete(pb);
     pb=NULL;
     return;
    }

 if (Draw)
    {
     Pf=EvalpF(*pb,1,BetterUbf);
     if ( Pf < (REAL) 0.1 )
        pb->Draw(XLim,WWidth,"Blue");
     if ( Pf >= (REAL)0.1 && Pf < (REAL) 0.2 )
        pb->Draw(XLim,WWidth,"Red");
     if ( Pf >= (REAL)0.2 && Pf < (REAL) 0.3 )
        pb->Draw(XLim,WWidth,"Green");
     if ( Pf >= (REAL)0.3 && Pf < (REAL) 0.4 )
        pb->Draw(XLim,WWidth,"Yellow");
     if ( Pf >= (REAL) 0.4 )
        pb->Draw(XLim,WWidth,"White");
    }

 // Termination criterion: the box and its objective range are not yet small enough.
 if ( NDiamGtEp (pb->X, Epsilon) && RDiam(pb->FX) > Epsilon )
    {
     pbltwork->Insert(pb, Inf(pb->FX));
    }
 else // the box goes to the final tree, after the remaining tests
    {
     if ( MTest )
        {
	 fgEvalH(f, pb->X, pb->FX, GX); NEvalG++;
	 if ( MonoTestFeas(pb, GX, NDelMT, XLim, WWidth, Draw) )
	    return;

	 // Try to improve the upper bound at the Baumann point, and use the
	 // centered form at that point to improve the bounds of the box.
	 pc = new BOX(*pb);
	 BaumannPoint(pb, GX, pc); 
	 fEvalH (f, pc->X, pc->FX); NEvalFx++;
	 if (Feasible(pc)==1 && UpdateBetterUbf (pb, BetterUbf, BetterPoint, f, NEvalFx, pc)) 
	   {
	    CutOffTest ( pbltwork, BetterUbf, NDelCTT, XLim, WWidth, Draw, "Brown" );
	    CutOffTest ( pbltfinal, BetterUbf, NDelCTT, XLim, WWidth, Draw, "Brown" );	  
	    cerr << NDelCTT << "\t" << pbltwork->Dimension() << "\t" << pbltfinal->Dimension() << endl;    
	   }
	 Intersection(pb->FX, pb->FX, (pc->FX + (GX)*(pb->X - pc->X))); 
	 delete pc;
	 pc = NULL;

	 if ( CutOffTheActualBox ( pb, BetterUbf, NDelMPT, XLim, WWidth, Draw)) 
	    return;
	}
     pbltfinal->Insert(pb, Inf(pb->FX));
    }
}

/*----------------------------------------------------------------------------*/
void DivideBox(const PBLTREE & pbltwork, const PBLTREE & pbltfinal, 
		PBOX & pb, INTERVAL_VECTOR & GX, PBOX & pc, 
               HTscalar_FctPtr f, const int & NDimDiv, 
               const int & NFrac, REAL & BetterUbf, VECTOR & BetterPoint,
               REAL Epsilon,
               const INTERVAL_VECTOR & XLim, INT WWidth, 
               bool Draw, bool MTest, bool EFatherGX,
	       int & NEvalFx, int & NEvalFX, int & NEvalG, 
	       int & NDelMPT, int & NDelMT, int & NDelCTT, int & NDelFathGX)
{
 int NSubBoxes;  	// number of subboxes
 PBOX pSubBox;
 PREAL Incr;            // width of the parts in each dimension
 PREAL WDim;		// widths of the dimensions
 PINT  IndWDim;         // indices of the dimensions ordered by width
 REAL  lb, ub;
 int i, j, Ind, Dim;

 Dim       = Dimension(*pb);
 NSubBoxes = (int) pow ((REAL)NFrac, (REAL)NDimDiv);
 
 Incr = new REAL[Dim];
 
 for ( i=0; i<Dim; i++ )
     Incr[i] = Diam( pb->X(i+1) ) /(REAL)NFrac;

 if (NDimDiv == Dim)    // divide all dimensions
    for ( i=0; i<NSubBoxes ;i++ )
	{
	 pSubBox = new BOX(true, *pb);
	 for ( j=0; j<Dim; j++ )
         {
          lb = Inf(pb->X(j+1)) +
               (REAL)((i/(int)pow((REAL)NFrac,(REAL)j))%NFrac)*Incr[j];
          ub = lb + Incr[j];
          pSubBox->X(j+1) = Hull(lb,ub);
         }    

	 CheckBox(pbltwork, pbltfinal, pSubBox, GX, pc, f, BetterUbf, BetterPoint, Epsilon,
	          XLim, WWidth, Draw, MTest, EFatherGX,
		  NEvalFx, NEvalFX, NEvalG,
	          NDelMPT, NDelMT, NDelCTT,  NDelFathGX);
      }
 else                   // divide the NDimDiv widest dimensions
    {
     VECTOR Widths(Dim);
     IndWDim = new int[Dim];
     CheckMemory(IndWDim);
     WDim    = new REAL[Dim];
     CheckMemory(WDim);

     Widths=Diam(pb->X);
     
     for (i=0;i<Dim;i++)
         {
          WDim[i]=Widths(i+1);
          IndWDim[i] = i;
         } 

     Bubble(IndWDim,WDim,Dim);

     for ( i=0; i<NSubBoxes ;i++ )
	 {
	  pSubBox = new BOX(true, *pb);
	  for ( j=0; j<NDimDiv; j++ )
          {
           Ind=IndWDim[j];
           lb = Inf(pb->X(Ind+1)) +
                (REAL)((i/(int)pow((REAL)NFrac,(REAL)j))%NFrac)*Incr[Ind];
           ub = lb + Incr[Ind];
           pSubBox->X(Ind+1) = Hull(lb,ub);
          }	

	  CheckBox(pbltwork, pbltfinal, pSubBox, GX, pc, f, BetterUbf, BetterPoint, Epsilon,
	           XLim, WWidth, Draw, MTest, EFatherGX,
		   NEvalFx, NEvalFX, NEvalG,
	           NDelMPT, NDelMT, NDelCTT, NDelFathGX);
       }	    
         
     delete [] WDim;
     delete [] IndWDim;
   }     

 delete [] Incr;
 if (Draw)
    pb->UnDraw();
 delete pb;
 pb=NULL;
}
