/*
 NestedLogit: interval branch-and-bound method for the continuous competitive
 single facility location and design problem with the nested logit choice
 rule.

 The new facility of the locating chain (chain 1) is described by its location
 (x1, x2) and its quality alpha. The program minimises the negative profit of
 the chain (see Problems.cpp) over the box given in the instance file, subject
 to the constraint that the facility is not too close to the demand points.

 The algorithm keeps a working tree and a final tree of boxes, both sorted by
 the lower bound of the objective. In each iteration the box with the smallest
 lower bound is taken from the working tree, the rejection tests are applied
 (cut-off, monotonicity and feasibility tests), the upper bound of the global
 minimum is updated, and the box is divided. Boxes that
 satisfy the termination criterion are moved to the final tree.

 Usage: NestedLogit -f <instance> [options]     (NestedLogit -h for the options)
 Output: a report on stderr and in <instance><options>.res, and the final
 boxes in <instance><options>_res.tclin.
*/

#include <math.h>
#include <time.h>
#include <values.h>
#include <iostream>
#include <fstream>
#include "I_util.h"
#include "Interval.h"
#include "Vector.h"
#include "IntervalVector.h"
#include "Matrix.h"
#include "Box.h"
#include "BLTree.h"
#include "Problems.h"
#include "ProcessBox.h"
#include "GORoutines.h"
#include "InputOutput.h"
#include <sys/times.h>
#include <unistd.h>

using namespace std;

/*----------------------------------------------------------------------------*/
/* Data of the problem, read from the instance file (see InputOutput.cpp).    */
/*----------------------------------------------------------------------------*/
int n, m, k;         // n - number of demand points, m - number of existing facilities,
                     // k - number of existing facilities of the locating chain
REAL b1, b2;         // parameters of the distance function
REAL a0, a1;         // parameters of the quality cost function
REAL c;              // income per unit of captured demand
MATRIX fi;           // fi(i,2): parameters of the cost of the land G1(x)
MATRIX p;            // locations of the demand points
VECTOR w;            // buying power at the demand points
VECTOR v;            // weight of the quality of the new facility at the demand points
MATRIX a, f;         // qualities and locations of the existing facilities
REAL ww;             // sum of the buying power
VECTOR t2, r2;       // precomputed terms of the probabilistic rule
VECTOR aJ1, b;       // aJ1: own chain term, b: competing chains term of the nested logit
VECTOR Minima;       // a point given in the instance (used by -ft)
REAL delta, lac, la1; // delta: forbidden radius parameter, lambda_c and lambda_1
REAL Beta = 1.0;     // distance decay: g(d) = exp(Beta*d)
bool fixvalues = 0;  // -fv

/* Colours of the hulls of the solution boxes in the graphical output. */
static const char* color[] = {
  "RosyBrown3", "IndianRed3", "IndianRed4", "sienna2", "burlywood3",
  "wheat2", "chocolate2", "chocolate4", "firebrick3", "firebrick4",
  "DodgerBlue3", "orange1", "DodgerBlue2", "DodgerBlue1", "DodgerBlue4",
  "SteelBlue1", "SteelBlue3", "SkyBlue2", "LightCyan3", "PaleTurquoise3",
  "CadetBlue3", "turquoise4", "SeaGreen3", "PaleGreen3", "SpringGreen3",
  "green4", "OliveDrab3", "DarkOliveGreen4", "khaki3", "khaki4",
  "yellow1", "yellow2", "yellow3", "gold1", "gold3",
  "goldenrod2", "RoyalBlue4", "blue3", "blue4", "bisque3",
  "LemonChiffon2", "honeydew3", "LavenderBlush3", "MistyRose3", "azure3",
  "SlateBlue3", "SlateBlue1", "DarkOrange2", "orchid3", "MediumPurple3",
  "RosyBrown3", "IndianRed3", "IndianRed4", "sienna2", "burlywood3",
  "wheat2", "chocolate2", "chocolate4", "firebrick3", "firebrick4",
  "orange1", "DodgerBlue1", "DodgerBlue2", "DodgerBlue3", "DodgerBlue4",
  "SteelBlue1", "SteelBlue3", "SkyBlue2", "LightCyan3", "PaleTurquoise3",
  "CadetBlue3", "turquoise4", "SeaGreen3", "PaleGreen3", "SpringGreen3" };
static const int NColors = sizeof(color)/sizeof(color[0]);

/*----------------------------------------------------------------------------*/
/* Write the demand points (with their forbidden circles) and the existing    */
/* facilities for the graphical output.                                        */
/*----------------------------------------------------------------------------*/
static void DrawDemandPoints(ostream & out, const INTERVAL_VECTOR & XLim, int WWidth)
{
 REAL FXDiam = Diam(XLim(1)), FYDiam = Diam(XLim(2));

 for (int i=1; i<=n; i++)
     {
      out << "PintaCircle" << endl;
      out << (p(i,1)-w(i)/(sqrt(b1)*delta)-Inf(XLim(1)))  / FXDiam * (REAL)(WWidth-1) << endl; 
      out << ((p(i,2)-w(i)/(sqrt(b2)*delta)-Inf(XLim(2))) / FYDiam * (REAL)(-1.0) + (REAL)1.0 ) *
             (REAL)(WWidth-1) << endl;
      out << (p(i,1)+w(i)/(sqrt(b1)*delta)-Inf(XLim(1)))  / FXDiam * (REAL)(WWidth-1) << endl; 
      out << ((p(i,2)+w(i)/(sqrt(b2)*delta)-Inf(XLim(2))) /FYDiam  * (REAL)(-1.0) + (REAL)1.0 ) *
             (REAL)(WWidth-1) << endl; 
      out << "Black" << endl;
      out << i << endl;
     }
}

static void DrawFacilities(ostream & out, const INTERVAL_VECTOR & XLim, int WWidth)
{
 REAL FXDiam = Diam(XLim(1)), FYDiam = Diam(XLim(2));

 for (int i=1; i<=m; i++)
     {
      out << "PintaCorner" << endl;
      out << (f(i,1)-Inf(XLim(1)))/FXDiam*(REAL)(WWidth-1) << endl; 
      out << ((f(i,2)-Inf(XLim(2)))/FYDiam*(REAL)(-1.0) + (REAL)1.0 ) *
             (REAL)(WWidth-1) << endl;
      if (i>k) out << "Red" << endl;          // competitor
      else out << "Darkgreen" << endl;        // own chain
      out << "2" << endl;
     }
}

/*----------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
 BLTREEP_VECTOR pbltv(2);      // pbltv(1): working tree, pbltv(2): final tree
 INTERVAL_VECTOR XLim, GX;     // search region, gradient inclusion
 VECTOR BetterPoint;           // best feasible point found
 PBOX pb=NULL;                 // the current box
 PBOX pc=NULL;                 // the box of the midpoint or Baumann point "c"
 REAL BetterUbf, Epsilon;      // upper bound of the global minimum, tolerance
 int NEvalFx, NEvalFX, NEvalG; // number of point, interval and gradient evaluations
 int NDelMPT, NDelCTT, NDelMT, NDelFathGX; // number of boxes deleted by the tests
 int NDim, NDimDiv, NFrac, WWidth, NBoxExt, i;
 bool MTest, FTest, Draw, EFatherGX, UpdBau;
 char NameProblem[256], fileres[256];
 struct tms tmsave1, tmsave2;
 double tp, ct = (double) sysconf(_SC_CLK_TCK);
 ofstream FResults, res;
 HTscalar_FctPtr ff=fFacLoc, ffD=fFacLocD;
 INTERVAL FB, FBvi;

 cerr.precision(16);
 cout.precision(16);

 GetParams (argc, argv, NameProblem, ff, NDim, XLim, Minima,
            Epsilon, NFrac, NDimDiv, MTest, EFatherGX, UpdBau,
            WWidth, FTest, Draw);

 Resize(GX,NDim);

 pb = new BOX(NDim);
 pc = new BOX(NDim);

 // -ft: evaluate the problem at the point Minima and stop.
 if (FTest)
    { 
     pb->X=Hull(Minima);
     fgEvalH(ffD,pb->X,pb->FX,GX); NEvalG++;
     if (Feasible(pb)==-1) cerr << "Infeasible.\n";
     cerr << *pb << endl;
     exit(0);
    }			

 // Header of the graphical output.
 if (Draw)
     {              
      cout << WWidth << endl;
      cout << WWidth << endl;
      cout << Inf(XLim(1)) << endl;
      cout << Sup(XLim(1)) << endl;
      cout << Inf(XLim(2)) << endl;
      cout << Sup(XLim(2)) << endl;
      cout << Epsilon << endl;
      cout << NameProblem << " problem" << endl;
      DrawDemandPoints(cout, XLim, WWidth);
     }

 times(&tmsave1);

 NEvalFx = NEvalFX = NEvalG = 0;
 NDelMPT = NDelCTT = NDelMT = NDelFathGX = 0;
 NBoxExt = 0;
 BetterUbf = DBL_MAX;
 Resize (BetterPoint,NDim);
 
 // Initial upper bound: the midpoint of the search region, if feasible.
 BetterPoint = Mid(XLim);
 pc->X = BetterPoint;
 fEvalH (ff, pc->X, pc->FX); NEvalFx++;
 if ( Feasible(pc) == 1)
   UpdateBetterUbf (pb, BetterUbf, BetterPoint, ff, NEvalFx, pc);
 else Clear(BetterPoint);
   
 if (!MTest) 
   InitBox ( pb, XLim, ff, NEvalFX );
 else
   InitBoxGrad ( pb, GX, XLim, ff, NEvalG );

 CheckBox ( pbltv(1), pbltv(2), pb, GX, pc, ff, BetterUbf, BetterPoint, Epsilon, 
            XLim, WWidth, Draw, MTest, EFatherGX,
            NEvalFx, NEvalFX, NEvalG, NDelMPT, NDelMT, NDelCTT, NDelFathGX);

 /*--------------------------------------------------------------------------*/
 /* Main loop: process the box with the smallest lower bound.                */
 /*--------------------------------------------------------------------------*/
 while ( !pbltv(1)->isEmpty() )
       {
        pbltv(1)->GetLwKeyBox(pb);	

	// A fresh copy, to have the right feasibility flags in pc.
	delete pc;
	pc = new BOX(*pb);
	 
        if (MTest)
	   { 
	    fgEvalH(ff, pb->X, pb->FX, GX); NEvalG++;

	    if ( MonoTestFeas(pb, GX, NDelMT, XLim, WWidth, Draw) )
	       {
        	NBoxExt++;
		continue;
	       }

	    if (UpdBau)
	      BaumannPoint(pb, GX, pc);  // pc is the Baumann point of pb
	    else
	      pc->X = Mid(pb->X);        // pc is the midpoint of pb

	    fEvalH (ff, pc->X, pc->FX); NEvalFx++;
	    Feasible(pc);
	    // The centered form: F(c) + G(X)*(X-c).
	    Intersection(pb->FX, pb->FX, (pc->FX + (GX)*(pb->X - pc->X))); 

	    if ( CutOffTheActualBox ( pb, BetterUbf, NDelMPT, XLim, WWidth, Draw)) 
	       {
        	NBoxExt++;
		continue;
	       }
	   }
        else
	   {
	    pc->X = Mid(pb->X);
	    if ( Feasible(pc) == 1 ) 
	      {
	       fEvalH (ff, pc->X, pc->FX); NEvalFx++;		   
	      }
	   }
	         
        // Update the upper bound and remove the boxes above it.
        if ( pc->feas==true )
	   if ( UpdateBetterUbf (pb, BetterUbf, BetterPoint, ff, NEvalFx, pc) )
             {
              CutOffTest ( pbltv(1), BetterUbf, NDelCTT, XLim, WWidth, Draw, "Brown" );
              CutOffTest ( pbltv(2), BetterUbf, NDelCTT, XLim, WWidth, Draw, "Brown" );
	      cerr << NDelCTT << "\t" << pbltv(1)->Dimension() << "\t" << pbltv(2)->Dimension() << endl;
             }

	DivideBox ( pbltv(1), pbltv(2), pb, GX, pc, ff, NDimDiv, NFrac, 
		       BetterUbf, BetterPoint, Epsilon, XLim, WWidth, Draw,
		       MTest, EFatherGX, NEvalFx, NEvalFX, NEvalG, NDelMPT, 
		       NDelMT, NDelCTT, NDelFathGX);
                   
	NBoxExt++; 

	// For large problems: report the time, and stop after 6 hours.
 	if (n>=1000 && NBoxExt>1000 && NBoxExt % 500 == 0)
	  {
	   times(&tmsave2); 
	   tp = tmsave2.tms_utime + tmsave2.tms_cutime -tmsave1.tms_utime - tmsave2.tms_cutime;
	   cerr << "Time: " << tp/ct << endl;
	  
	   if (tp/ct > 6*3600) 
	     { 
	      cerr << "Final time: " << tp/ct << endl;
	      cerr << "Upper bound: " << BetterUbf << endl;
	      cerr << "Lower bound: " << pbltv(1)->MinlbF() << " or " << pbltv(2)->MinlbF() << endl; 
	      exit(1);
             }
	  }
       } 

 delete pc;
 pc = NULL;

 times(&tmsave2); 
 tp = tmsave2.tms_utime + tmsave2.tms_cutime -tmsave1.tms_utime - tmsave2.tms_cutime;

 /*--------------------------------------------------------------------------*/
 /* Reports.                                                                 */
 /*--------------------------------------------------------------------------*/
 INTERVAL_MATRIX Results;      // hulls of the groups of final boxes
 Resize(Results,10,NDim);
 Clear(Results);

 OutPutStat (cerr, pbltv, Epsilon, NFrac, NDimDiv, BetterPoint, BetterUbf,
             tp/ct, NEvalFx, NEvalFX, NEvalG, NameProblem, MTest,
             NDelMPT, NDelCTT, NDelMT, NDelFathGX);

 // The name of the result files: the instance name followed by the options.
 strcpy(fileres,NameProblem);
 for (i=3;i<argc;i++)
     strcat(fileres,argv[i]);
 strcat(fileres,".res");
 res.open(fileres, ios::out);
 res.precision(16);

 OutPutStat (res, pbltv, Epsilon, NFrac, NDimDiv, BetterPoint, BetterUbf,
             tp/ct, NEvalFx, NEvalFX, NEvalG, NameProblem, MTest,
             NDelMPT, NDelCTT, NDelMT, NDelFathGX);
 
 cerr << "Final time: " << tp/ct << endl;
 res << "Final time: " << tp/ct << endl;

 // The objective at the best point with both choice rules.
 fEvalH(&fFacLocD, Hull(BetterPoint), FB); 
 fEvalH(&fFacLocViejoD, Hull(BetterPoint), FBvi);
 cerr << "At the best point using Nested-Logit: " << Inf(FB) << endl;
 cerr << "At the best point using Probablist: " << Inf(FBvi) << endl;
 res << "At the best point using Nested-Logit: " << Inf(FB) << endl;
 res << "At the best point using Probablist: " << Inf(FBvi) << endl;

 pbltv(2) = ComputeResults(pbltv(2),Results);
 cerr << "Hull: " << endl;
 res << "Hull: " << endl;
 i=1;
 while ( Inf( Results(i,3) ) > 0.0 )
     {
      cerr << "([" << Inf(Results(i,1)) << "," << Sup(Results(i,1)) << "] ; [" << Inf(Results(i,2)) << "," << Sup(Results(i,2)) << "] ; [" << Inf(Results(i,3)) << "," << Sup(Results(i,3)) << "])" << endl;
      res << "([" << Inf(Results(i,1)) << "," << Sup(Results(i,1)) << "] ; [" << Inf(Results(i,2)) << "," << Sup(Results(i,2)) << "] ; [" << Inf(Results(i,3)) << "," << Sup(Results(i,3)) << "])" << endl;
      i++;
     }
 res.close();
 
 // The hulls of the solution boxes, the demand points and the facilities
 // for the graphical output.
 if ( !(pbltv(2)->isEmpty()) )
     {              
      strcpy(fileres,NameProblem);
      for (i=3;i<argc;i++)
	   strcat(fileres,argv[i]);
      strcat(fileres,"_res.tclin");
      cout << fileres << endl;
      res.open(fileres, ios::out);
      res.precision(16);
      res << WWidth << endl;
      res << WWidth << endl;
      res << Inf(XLim(1)) << endl;
      res << Sup(XLim(1)) << endl;
      res << Inf(XLim(2)) << endl;
      res << Sup(XLim(2)) << endl;
      res << Epsilon << endl;
      res << NameProblem  << endl;

      pc = new BOX(NDim); i=1;
      while ( Inf( Results(i,3) ) > 0.0 )
         {
	  pc->X = Row(Results,i++);
	  pc->Draw(XLim,WWidth,(PCHAR)color[(2*i) % NColors], res);
	 }
      DrawDemandPoints(res, XLim, WWidth);
      DrawFacilities(res, XLim, WWidth);
     }

 if (ExistArg("-FEst",argc,argv))           
    {
     FResults.open(GetArg("-FEst",argc,argv), ios::app);
     OutPutStat (FResults, pbltv, Epsilon, NFrac, NDimDiv, BetterPoint, BetterUbf,
                 tp/ct, NEvalFx, NEvalFX, NEvalG, NameProblem, MTest,
                 NDelMPT, NDelCTT, NDelMT, NDelFathGX);
     FResults.close();
    }
 }
