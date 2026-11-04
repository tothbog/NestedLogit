#include <math.h>
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h> 
#include <stdlib.h>

typedef char * PCHAR; 
typedef PCHAR * PPCHAR; 

using namespace std;

/*----------------------------------------------------------------------------*/
/* Look for the string 'str_to_find' in 'arg'.                                */
/* If it is found, return the next string in arg.                             */
/* In other case return NULL.*/
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
/* Return 1 if the strind 'str_to_find' is in the command line 'arg'.*/
/*----------------------------------------------------------------------------*/
bool ExistArg (PCHAR str_to_find, int narg, PPCHAR arg) 
{
 int i;

 for ( i=0; i<narg; i++)
     if ( strstr(arg[i], str_to_find) && strstr(str_to_find, arg[i]) ) 
        return true;
 return false;
}
/*----------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{

  char foo[256], ProblemName[256];
  string fname("r-");
  int num;          // number of problems
  int n, m, k;       // n - number of demand points,  m - number of existing facilities,
  	             // k - lenght of the chain of existing facilities
  MATRIX bounds(3,2); // bounds for x and alpha
  double b1, b2, pp; // parameters for the distance function
  double l;          // parameter for the g functions
  double a0, a1;     // parameters for the G(alpha) functions
  double la1, lac;     // parameters for the Nested Logit 
  double c;	     // parameter for F
  int parts;
  
  srand48(time(0));  // inicialize the random generator with the actual time
  
  //n = 10 + trunc(40 * drand48());
  //m = 3 + trunc(7 * drand48());
  //k = trunc(m/2 * drand48());
  
  if ( argc <4 )
    {
     cerr << "Usage: ProblemGenerator NumOfProblemToGenerate n m k [UpperBoundOfTheSearchingRegion]\n\n";
     exit(1);
    }
 
  if (!ExistArg("-p",argc,argv))
    parts = 1;
 else
    parts = (REAL)atoi(GetArg("-p",argc,argv));
 
  num = atoi(argv[1]);
  n = atoi(argv[2]);
  m = atoi(argv[3]);
  k = atoi(argv[4]);
  fname += argv[2];
  fname += "-";
  fname += argv[3];
  fname += "-";
  fname += argv[4];
  fname += "-";
  
  bounds(1, 1) = bounds(2, 1) = 0.0; bounds(3, 1) = 0.5;
  if (argc > 5)
     bounds(1,2) = bounds(2,2) = bounds(3, 2) = atoi(argv[5]);
  else   
      bounds(1, 2) = bounds(2, 2) = 10.0;
  bounds(3, 2) = 5;
      
  //cerr << fname << endl; exit(1);
  
  VECTOR chain(k);
  MATRIX fi(n, 2);    // parameters for the Psi functions 
  MATRIX f(m, 2);     // f - the location of the existing facilities
  MATRIX p(n, 2);     // p - the location of the demand points
  VECTOR w(n);        // w - buying power at the demand points
  VECTOR v(n);        // v - weight of attractiveness to x of the demand points
  MATRIX a(1, m);     // attractiveness of demand points to the existing facilities
  VECTOR Minima(3);
  double Min; 
  int i, j;
  

  for (int lll=0; lll<num; lll++)
  {
      //for (int p=1; p<=parts; p++)
      sprintf(foo, "%d-11", lll);
      ofstream fout((fname + foo).c_str());
      if (!fout.good()) {cerr << "Could not open file: " << foo << endl; exit(1);}
      strcpy(ProblemName, (fname + foo).c_str());
      
      b1 = 1;// + drand48();
      b2 = 1;// + drand48();
      la1 = 1;// + drand48();
      lac = 1;// + drand48();
      pp = l = 2;
      a0 = 7 + 2 * drand48();
      a1 = 6 + 0.5 * drand48();
      c = 10 + 2 + 2*drand48(); // large problems are without the 2+2*, but for table 4 I had to change it.

      fout << "Name_of_problem			" << ProblemName << endl;
      fout << "num_of_demand_points(n)		" << n << endl;
      fout << "num_of_existing_facilities(m)	" << m << endl;
      fout << "number_of_owned_fac(k)		" << k << endl;
      fout << "Bounds_for_the_variables_x1_x2_alpha:\n";
      fout << bounds(1, 1) << "\t" << bounds(1, 2) << "\n";
      fout << bounds(2, 1) << "\t" << bounds(2, 2) << "\n";
      fout << bounds(3, 1) << "\t" << bounds(3, 2) << "\n\n";
      fout << "Infeasible_radius_parameter(delta): " << (3.0*sqrt(n+50.0)-20.0) << endl << endl;
      fout << "Parameters_for_distance_function_(b1*(x1-x2)^p+b2*(x2-y2)^p)^1/p:" << endl;
      fout << "b1	" << b1 << endl;
      fout << "b2	" << b2 << endl << endl;
//      fout /* << "p	" */ << pp << endl << endl;
//      fout /* << "Parameter_for_g_function(lambda__or__beta__in_hodson): " */ << l << endl << endl;
      fout << "Parameters_for_Nested_Logit(lambda1,lambda_c):" << endl;
      fout << "a0	" << a0 << endl;
      fout << "a1	" << a1 << endl << endl;
      fout << "Parameters_for_G(alpha)_functions:" << endl;
      fout << "la1 " << la1 << endl;
      fout << "lac  " << lac << endl << endl;


      fout << "Parameter_for_F: "   << c << endl << endl;  


      fout << "Parameters_for_Psi_functions:\n";
      for (i=1; i <=n; i++)
      {
	  fi(i, 1) = 2;
	  fi(i, 2) = 0.5 + 1.0 * drand48(); // 0.5 + 1.5 * drand48(); was

	  fout << fi(i, 1) << "\t" << fi(i, 2) << endl;
      }
      fout << endl << endl << "location_of_the_existing_facilities(f[m][2]):\n"; 
      for (j=1; j <=m; j++)
      {
	  f(j, 1) = bounds(1,1) + bounds(1,2) * drand48();
	  f(j, 2) = bounds(2,1) + bounds(2,2) * drand48(); 

	  fout << f(j, 1) << "\t" << f(j, 2) << endl;
      }      
      fout << endl << endl << "location_of_the_demand_points(p[n][2]):\n";
      for (i=1; i <=n; i++)
      {
	  p(i, 1) = bounds(1,1) + bounds(1,2) * drand48();
	  p(i, 2) = bounds(2,1) + bounds(2,2) * drand48(); 

	  fout << p(i, 1) << "\t" << p(i, 2) << endl;
      }
      fout << endl << endl << "buying_power_at_the_demand_points(w[n]):\n"; 
      for (i=1; i <=n; i++) // !!!!! NEW STYLE!!! max 10 for 100 demand points
      {
	  w(i) = 100/sqrt(n) * drand48();
	  fout << w(i) << endl;
      }
      fout << "\n\n" << "weight_of_attractiveness_to_x_of_the_demand_points(v[n]--gamma):\n"; 
      for (i=1; i <=n; i++)
      {
	  v(i) = 0.75 + 0.5 * drand48();
	  fout << v(i) << endl;
      }

      fout << endl << endl << "attractiveness_of_demand_points_to_the_existing_facilities(a[m]--alfa):\n";
      for (i=1; i <=1; i++) // !!!!! NEW STYLE!!! ONLY ONE LINE IS WRITTEN!
      {
	  for (j=1; j <=m; j++)
	  {
	      if (i>1)
	      	a(i, j) = a(1,j);
	      else
	      	a(i, j) = bounds(3,1) + bounds(3,2) * drand48();
	      fout << a(i, j) << "\t"; 
	  }   
	  fout << endl;   
      }
      fout << endl;   



      Minima(1) = bounds(1,1) + bounds(1,2) * drand48();
      Minima(2) = bounds(2,1) + bounds(2,2) * drand48();
      Minima(3) = bounds(3,1) + bounds(3,2) * drand48();
      Min = 100.0;

       fout << "Minima:\n" << Minima(1) << endl << Minima(2) << endl << Minima(3) << endl << endl << endl; 
       fout << "Minimo:\n" << Min;
      
      fout.close();
  } 
}
