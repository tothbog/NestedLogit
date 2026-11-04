#include <stdlib.h>
#include <fstream>
#include <string>
#include "InputOutput.h"

extern  int n, m, k;       // n - number of demand points, m - number of existing facilities,
  	                   // k - number of existing facilities of the locating chain (the first k)
extern  REAL b1, b2;       // parameters of the distance function
extern  REAL a0, a1;       // parameters of the quality cost function
extern  REAL c;            // income per unit of captured demand
extern  REAL ww;           // sum of the buying power
extern  VECTOR w;	   // buying power at the demand points
extern  VECTOR v;          // weight of the quality of the new facility at the demand points
extern  MATRIX a;          // a(1,j): quality of existing facility j
extern  MATRIX f;          // locations of the existing facilities
extern  MATRIX fi;         // fi(i,2): parameters of the cost of the land G1(x)
extern  MATRIX p;          // locations of the demand points
extern  VECTOR t2,r2;      // probabilistic rule: w(i)*attraction of the competitors, attraction of all
extern  VECTOR aJ1, b;     // nested logit: own chain term, aggregated competitor term
extern  REAL delta, lac, la1; // delta: forbidden radius parameter, lambda_c and lambda_1 of the nested logit
extern  REAL Beta;         // distance decay parameter
extern  bool fixvalues;    // -fv: fix the land cost parameters and the weights v

/*----------------------------------------------------------------------------*/
void InputError()
{
 cerr << "Options are:" << endl;
 cerr << "\t[ -h          : To show this help ]" << endl;
 cerr << "\t  -f   <file> : Name of the file containing the problem" << endl;
 cerr << "\t[ -ep  <real> : Epsilon, tolerance of the termination rule (default 1e-3) ]" << endl;
 cerr << "\t[ -mt         : Monotonicity test ]" << endl;
 cerr << "\t[ -fv         : Fix fi(i,2)=1.25 and v(i)=1 for all demand points ]" << endl;
 cerr << "\t[ -ndv  <int> : Number of parts per divided dimension (default 2) ]" << endl;
 cerr << "\t[ -ndi  <int> : Number of dimensions to divide, 0 for all (default 2) ]" << endl;
 cerr << "\t[ -ef         : Bound the subboxes with the father's gradient (only with -mt) ]" << endl;
 cerr << "\t[ -ub         : Update the upper bound at the Baumann point (only with -mt) ]" << endl;
 cerr << "\t[ -vi         : Use the probabilistic choice rule instead of the nested logit ]" << endl;
 cerr << "\t[ -ft         : Print the input data and test the function at Minima ]" << endl;
 cerr << "\t[ -FEst <file>: Append the statistics to file ]" << endl;
 cerr << "\t[ -w    <int> : Width of the graphical output window (default 800) ]" << endl;
 cerr << "\t[ -tcl        : Graphical output to stdout ]" << endl;
 exit(0);
}                  

/*----------------------------------------------------------------------------*/
/* Two boxes are near if the sum of their distances along the coordinates is  */
/* not larger than the sum of their widths.                                    */
/*----------------------------------------------------------------------------*/
static INT NearBy(INTERVAL_VECTOR & X,INTERVAL_VECTOR & Y)
{
 int Dim = Dimension(X);
 REAL dist = 0.0, diam = 0.0;
 for(int i=1; i<=Dim; i++)
    {
     dist += Distance(X(i),Y(i));
     diam += Diam(X(i))+Diam(Y(i));
    }
 
 if (dist <= diam) 
    return 1;
 else
    return 0;
}

/*----------------------------------------------------------------------------*/
/* Merge the rows of Results (hulls) that became near to each other.           */
/*----------------------------------------------------------------------------*/
static void Pack(INTERVAL_MATRIX & Results, INT & sols)
{
 int Dim = ColDimension(Results);
 INTERVAL_VECTOR Y(Dim), X(Dim);
 
  for (int j=1; j<=sols;j++)
     {
      X = Row(Results,j);
      for(int i=j+1; i<=sols; i++)
	 {
	  Y = Row(Results,i);
	  if ( NearBy(X,Y) )
	    {
	     Y = Hull(X,Y);
	     X = Y;
	     SetRow(Results,j,Y);
	     for (int k=i; k<sols;k++)
		 SetRow(Results,k,Row(Results,k+1));
	     Clear(Y);
	     SetRow(Results,sols,Y);
	     sols--;i--;
	     continue;
	    }
	 }
     }
}

/*----------------------------------------------------------------------------*/
/* Empty the final tree pbltv, and collect the hulls of the groups of nearby   */
/* boxes in the rows of Results. The boxes are moved to a new tree, which is   */
/* returned.                                                                   */
/*----------------------------------------------------------------------------*/
PBLTREE ComputeResults(PBLTREE pbltv, INTERVAL_MATRIX & Results)
{
 int Dim = ColDimension(Results);
 PBOX pb = NULL;
 INTERVAL_VECTOR Y(Dim);
 INTERVAL_MATRIX T;
 int i,sols=0;
 bool neww;
 PBLTREE temp=NULL;
 temp = new BLTREE;

 Clear(Results);
 while ( !pbltv->isEmpty() )
     {
      pbltv->GetLwKeyBox(pb);
      neww=true;
      for(i=1; i<=sols; i++)
	 {
	  Y = Row(Results,i);
	  if ( NearBy(pb->X,Y) )
	    {
	     Y = Hull(pb->X,Y);
	     SetRow(Results,i,Y);
	     neww=false;
	    }
	 }
      if (sols > 1) Pack(Results, sols); 
      if (neww) 
         if (sols < RowDimension(Results)) 
           SetRow(Results,++sols,pb->X);
         else  // enlarge Results
	   {
	    Resize(T,sols,Dim); T=Results;
	    Resize(Results,sols+10,Dim);
	    for(i=1; i<=sols; i++)
	      SetRow(Results,i,Row(T,i));
	    SetRow(Results,++sols,pb->X);
	   }

      temp->Insert(pb, Inf(pb->FX));
     }
 return temp;
}

/*----------------------------------------------------------------------------*/
VOID OutPutStat (ostream & FileOut, BLTREEP_VECTOR & pbltv,
                 REAL Epsilon, int NFrac, int NDimDiv,
                 const VECTOR & BetterPoint, REAL BetterUbf,
                 double time, int NEvalFx, int NEvalFX, int NEvalG,
                 PCHAR NameProblem, bool MTest,
                 int NDelMPT, int NDelCTT, int NDelMT, int NDelFathGX)
{
 int Dim = Dimension(BetterPoint);
 int NFinalBoxes, MaxBoxesWTree;
 
 FileOut << "-----------------------------------------------------------"
         << endl;
 FileOut << "---                    NestedLogit-0.0.3                ---"
         << endl;
 FileOut << "-----------------------------------------------------------"
         << endl;
 FileOut << "Function                              = " << NameProblem << endl;
 FileOut << "Epsilon                               = " << Epsilon     << endl;
 FileOut << "Number of Divisions per Dimension     = " << NFrac-1     << endl;
 FileOut << "Number of divided dimensions one step = " << NDimDiv     << endl;
 FileOut << "Monocity Test                         = " << MTest       << endl;
 FileOut << "Boxes deleted by MidPoint     Test    = " << NDelMPT	  << endl;
 FileOut << "Boxes deleted by CutOff       Test    = " << NDelCTT	  << endl;
 FileOut << "Boxes deleted by Monotonicity Test    = " << NDelMT	  << endl; 
 FileOut << "Boxes deleted by EvalFatherGX Test    = " << NDelFathGX	  << endl;
 FileOut << "BetterPoint                           = " << BetterPoint << endl;
 FileOut << "f* in                                 [ " << pbltv(2)->MinlbF() << ", " << BetterUbf 
         << " ]" << endl;
 FileOut << "Execution time                        = " << time << endl; 
 FileOut << "Number of point    evaluations        = " << NEvalFx     << endl;
 FileOut << "Number of interval evaluations        = " << NEvalFX     << endl;
 FileOut << "Number of gradient evaluations        = " << NEvalG      << endl;
 FileOut << "Evaluation index (EvalF + n * EvalF') = " << NEvalFx+NEvalFX+ Dim*NEvalG << endl;

 NFinalBoxes   = pbltv(2)->Dimension();
 MaxBoxesWTree = pbltv(1)->MaxDimension();
 FileOut << "Number of boxes in the final tree     = " << NFinalBoxes << endl;
 FileOut << "Max number of nodes in work tree      = "
         << MaxBoxesWTree << endl;
    
 // One line for a LaTeX table:
 // time & point evals & interval evals & gradient evals & evaluation index &
 // final boxes & max work tree size & 0
 char S[256];
 if (MTest)   
     {
      sprintf(S,"%4.2f & %8d & %8d & %8d & %8d & %8d & %8d & %8d \\\\\n",
                  time, 
                  NEvalFx, NEvalFX, NEvalG, NEvalFx+NEvalFX+ Dim*NEvalG, NFinalBoxes,
		          MaxBoxesWTree, 0); 
      FileOut << S;
     }
 else
     {
       sprintf(S,"%4.2f & %8d & %8d & %8d & %8d & %8d \\\\\n",
                  time, NEvalFx, NEvalFX,
		          NFinalBoxes, MaxBoxesWTree, 0); 
      FileOut << S;
     }
		     
 FileOut << "-----------------------------------------------------------"
         << endl;
}

/*---------------------------------------------------------------------------*/
/* Read the instance file (see the README for its format) and precompute the */
/* terms of the objective that do not depend on the new facility.            */
/*---------------------------------------------------------------------------*/
void ReadProblem (PCHAR FileName, PCHAR NameProblem, INTERVAL_VECTOR & XLim, 
                  VECTOR & Minima, bool FTest)
{
 ifstream FIn;
 FIn.open(FileName, ios::in );
 if ( !FIn.good() )
    {
     cerr << "Unable to open file " << FileName << endl;
     exit(1);
    }

  char foo[256];
  int i,j;
  REAL Min;

  Resize(XLim, 3);

  // Each value is preceded by a label, which is skipped (read into foo).
  FIn >> foo >> NameProblem; strcpy(NameProblem, FileName);
  FIn >> foo >> n;
  FIn >> foo >> m;
  FIn >> foo >> k; 
  FIn >> foo >> XLim;
  FIn >> foo >> delta;
  FIn >> foo;
  FIn >> foo >> b1;
  FIn >> foo >> b2;
  FIn >> foo;
  FIn >> foo >> a0;
  FIn >> foo >> a1;
  FIn >> foo;
  FIn >> foo >> la1;
  FIn >> foo >> lac;
  FIn >> foo >> c;

  Resize(fi, n, 2);
  Resize(f, m, 2);
  Resize(p, n, 2);
  Resize(w,n);
  Resize(v, n);
  Resize(a, 1, m);
  Resize(t2,n);
  Resize(r2,n);
  Resize(aJ1,n);
  Resize(b,n);  

  FIn >> foo >> fi;
  FIn >> foo >> f; 
  FIn >> foo >> p;
  FIn >> foo >> w; 
  FIn >> foo >> v; 
  FIn >> foo >> a;

  if (fixvalues) {
    for(i=1;i<=n;i++) {
        fi(i,2) = 1.25;
        v(i) = 1;
    }
  }

  // A point (used only by -ft) and a value, kept for compatibility of the format.
  Resize(Minima,3);
  FIn >> foo >> Minima;
  FIn >> foo >> Min;

  if (FTest) {  
  cerr << "Name of the problem: " << NameProblem << endl;
  cerr << "n = " << n << endl;
  cerr << "m = " << m << endl;
  cerr << "k = " << k << endl;  
  cerr << "bounds =\n" << XLim << endl << endl;
  cerr << "delta= " << delta << endl;
  cerr << "b1 = " << b1 << endl;
  cerr << "b2 = " << b2 << endl;
  cerr << "a0 = " << a0 << endl;
  cerr << "a1 = " << a1 << endl;
  cerr << "la1 = " << la1 << endl;
  cerr << "lac = " << lac << endl;
  cerr << "c = " << c << endl;
  cerr << "f =\n" << f << endl << endl;
  cerr << "p =\n" << p << endl << endl;
  cerr << "w =\n" << w << endl << endl;
  cerr << "v =\n" << v << endl << endl;
  cerr << "a =\n" << a << endl << endl; 
 }

  // Attraction a(1,j)/g(d_ij) of the existing facilities, g(d) = exp(Beta*d).
  // Facilities 1..k belong to the locating chain, k+1..m to the competitors.
  REAL gdij;
  ww = 0.0; 
  for (i=1; i<=n; i++)
    {
      t2(i) = 0; r2(i) = 0; aJ1(i) = 0; b(i) = 0;
      for (j=1; j<=m; j++)
	{
	  gdij = exp(Beta*sqrt(b1 * pow(f(j,1)-p(i,1), 2) + b2 * pow(f(j,2)-p(i,2), 2)) + 1e-255); 
	  if (j > k) {
	    t2(i) += a(1,j)/gdij;                  // competitors
	    b(i)  += pow(a(1,j)/gdij,lac);
	  } else {
	    aJ1(i) += pow(a(1,j)/gdij,la1);        // own chain
	  }
	  r2(i) += a(1,j)/gdij;	                   // all existing facilities
	} 
      t2(i) *= w(i);
      ww += w(i);

      b(i) = pow(b(i),1.0/lac);
    }

  if (FTest) {
     cerr << "aJ1 =\n"  << aJ1 << endl << endl;
     cerr << "b =\n" << b << endl << endl;
     cerr << "t2 =\n"  << t2 << endl << endl;
     cerr << "r2 =\n" << r2 << endl << endl;
     cerr << "ww = " << ww << endl;
   }

 FIn.close();
}

/*----------------------------------------------------------------------------*/
/* Look for the string 'str_to_find' in 'arg'.                                */
/* If it is found, return the next string in arg, otherwise return NULL.      */
/*----------------------------------------------------------------------------*/
PCHAR GetArg (PCHAR str_to_find, int narg, PPCHAR arg) 
{
 int i;                              

 for ( i=0; i<narg; i++ )                
     if ( strstr(arg[i], str_to_find) && strstr(str_to_find, arg[i]) ) 
        return arg[i+1];          
 return NULL;           
}

/*----------------------------------------------------------------------------*/
/* Return true if the string 'str_to_find' is in the command line 'arg'.      */
/*----------------------------------------------------------------------------*/
bool ExistArg (PCHAR str_to_find, int narg, PPCHAR arg) 
{
 int i;

 for ( i=0; i<narg; i++)
     if ( strstr(arg[i], str_to_find) && strstr(str_to_find, arg[i]) ) 
        return true;
 return false;
}

/*---------------------------------------------------------------------------*/
void GetParams (int argc, PPCHAR argv, PCHAR NameProblem, HTscalar_FctPtr & ff, 
                int & NDim, INTERVAL_VECTOR & XLim, VECTOR & Minima,
                REAL & Epsilon, int & NFrac, int & NDimDiv,
                bool & MTest, bool & EFatherGX, bool & UpdBau,
		int & WWidth, bool & FTest, bool & Draw)
{
 PCHAR FileProblem;
 
 if (ExistArg("-h",argc,argv))
    InputError();   

 // Problem file
 if (!ExistArg("-f",argc,argv))
    {
     cerr << "Please specify the problem." << endl;
     InputError();
    }
 else
    FileProblem = GetArg("-f",argc,argv);     
 
 // Probabilistic choice rule instead of the nested logit
 if (ExistArg("-vi",argc,argv))
     ff = fFacLocViejo;

 // Test the function at Minima
 FTest = ExistArg("-ft",argc,argv);

 // Fix the land cost parameters fi(i,2) and the weights v(i)
 fixvalues = ExistArg("-fv",argc,argv);
 
 ReadProblem (FileProblem, NameProblem, XLim, Minima, FTest);
 NDim=Dimension(XLim);

 // Tolerance of the termination rule
 if (!ExistArg("-ep",argc,argv))
    Epsilon = 1e-3;
 else
    Epsilon = (REAL)atof(GetArg("-ep",argc,argv));

 // Number of parts in each divided dimension
 if (ExistArg("-ndv",argc,argv))
    {
     NFrac = atoi(GetArg("-ndv",argc,argv));      
     if ( NFrac <= 0)
        {
         cerr << "-ndv must be greater than 0" << endl;
         exit(1);
	}
    }	     
 else
    NFrac=2;  

 // Number of dimensions to divide    
 if (ExistArg("-ndi",argc,argv))
    {
     NDimDiv = atoi(GetArg("-ndi",argc,argv));
     if ( NDimDiv > NDim || NDimDiv <= 0 )
        NDimDiv=NDim;
    }
 else
    NDimDiv=2;

 // Monotonicity test
 MTest = ExistArg("-mt",argc,argv);

 // Bound the subboxes using the gradient of the father box
 EFatherGX = ExistArg("-ef",argc,argv);
 if (EFatherGX && !MTest)
    {
     cerr << "With -ef, -mt must be used.\n";
     exit(1);
    } 

 // Update the upper bound at the Baumann point
 UpdBau = ExistArg("-ub",argc,argv);
 if (UpdBau && !MTest)
    {
     cerr << "With -ub, -mt must be used.\n";
     exit(1);
    }

 // Width of the graphical output window
 if (ExistArg("-w",argc,argv))
     WWidth = atoi(GetArg("-w",argc,argv));
 else
     WWidth = 800;

 // Graphical output
 Draw = ExistArg("-tcl",argc,argv);
}
