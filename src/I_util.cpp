//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//                      Copyright (c) 1994 Rolf Hammer
//
// This program/module is free software for non-commercial use. For details
// on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
// This program/module is distributed WITHOUT ANY WARRANTY. For details,
// see the "Disclaimer / Legal Matters" of the book (page iv).
//
//============================================================================
// File: I_util (implementation)
// Purpose: Utilities of types 'INTERVAL', 'INTERVAL_VECTOR'
// Global functions:
//    in()      : contained-in relation for a real,
//                contained-in-the-interior relation for two intervals,
//                contained-in-the-interior relation for two interval
//                vectors and for an integer in an interval vector
//    Blow()    : epsilon inflation for an inteval,
//                epsilon inflation for an interval vector
//    Disjoint(): test for disjointedness of two intervals,
//                test for disjointedness of two interval vectors
//    AbsMin()  : smallest absolute value of an interval
//    AbsMax()  : greatest absolute value of an interval
//    RDiam()   : relative diameter of an interval (not the same as RelDiam!)
//    MaxRDiam(): to get the maximum of the relative diameters of the
//                components of an interval vector.
//    UlpAcc()  : to check whether the width of an interval is less than a
//                certain number of ulps (ulp = unit in the last place of
//                the mantissa),
//                to check whether all components of an interval vector
//                have width less than a certain number of ulps
//----------------------------------------------------------------------------
#include "Interval.h"	// Interval mathematic functions
#include "I_util.h"

/*----------------------------------------------------------------------------*/
bool in ( REAL x, const INTERVAL & y )                       
{ 
 return ( (Inf(y) <= x) && (x <= Sup(y)) ); 
}     

/*----------------------------------------------------------------------------*/
bool in ( const INTERVAL & x, const INTERVAL & y )      
{                                        
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

/*----------------------------------------------------------------------------*/
bool in ( INTERVAL_VECTOR & x, INTERVAL_VECTOR& y )
{                                          
 INT i = 1, n = Dimension(x); 
 bool ok = true;

 if ( n != Dimension(y) )
    return false;

 while (ok && i <= n) 
       { 
        ok = in( x(i),y(i) ); 
        i++; 
       }
 return ok;
}

/*----------------------------------------------------------------------------*/
bool in ( INT x, INTERVAL_VECTOR & y )     
{                                
  INT  i = 1, n = Dimension(y);
  bool ok = true;
  REAL d = x;

  while (ok && i <= n) 
        { 
         ok = in( d,y(i) ); 
         i++; 
        }
  return ok;
}

/*----------------------------------------------------------------------------*/
bool in ( const VECTOR & x, const INTERVAL_VECTOR & y )  
{                                
  INT  i = 1, n = Dimension(x); 
  bool ok = true;

  if ( n != Dimension(y) )
     return false;

  while (ok && i <= n) 
        { 
         ok = in( x(i),y(i) ); 
         i++;
        }
  return ok;
}

/*----------------------------------------------------------------------------*/
/*void rnd ( dotprecision& d, INTERVAL& x )    // Round dotprecision to interval
{                                            //-------------------------------
  real Lower, Upper;

  rnd(d,Lower,Upper);         // Rounding to interval bounds
  x = _INTERVAL(Lower,Upper);
}
*/

/*----------------------------------------------------------------------------*/
INTERVAL Blow ( INTERVAL x, REAL eps )                    // Epsilon inflation
{                                                         //------------------
  x =  (1.0+eps)*x - eps*x;
  return INTERVAL(Pred(Inf(x)),Succ(Sup(x)));
}

/*----------------------------------------------------------------------------*/
INTERVAL_VECTOR Blow ( INTERVAL_VECTOR& x, REAL eps )                     // Epsilon deflation
{                                                         //------------------
  INT             i, n = Dimension(x);
  INTERVAL_VECTOR h(n);
  
  for (i = 1; i <= n; i++) h(i) = Blow(x(i),eps);
  return h;
}

/*----------------------------------------------------------------------------*/
INT Disjoint ( INTERVAL a, INTERVAL b )             // Test for disjointedness
{                                                   //------------------------
  return (Inf(a) > Sup(b)) || (Inf(b) > Sup(a));
}

/*----------------------------------------------------------------------------*/
INT Disjoint ( INTERVAL_VECTOR& a, INTERVAL_VECTOR& b )
                                                    // Test for disjointedness
{                                                   //------------------------
  INT al = 1, au = Dimension(a);
  INT disjointed, i;

  i = al;
  disjointed = 0;
  do {
    if (Disjoint( a(i),b(i) ))
      disjointed = 1;
    else
      i++;
  } while ( !(disjointed || i > au) );
  return disjointed;
}

/*----------------------------------------------------------------------------*/
REAL AbsMin ( INTERVAL x )                              // Absolute minimum of
{                                                       // an interval
  if ( in(0.0,x) )                                      //--------------------
    return 0.0;
  else if (Inf(x) > 0.0)
    return Inf(x);
  else
    return -Sup(x);
}

/*----------------------------------------------------------------------------*/
REAL AbsMax ( INTERVAL x )                              // Absolute maximum of
{                                                       // an interval
//  REAL a, b;                                            //--------------------
 
//  a = fabs(Inf(x));
//  b = fabs(Sup(x));

  if (Inf(x) >= 0) return Sup(x);
  else if (Sup(x) <= 0) return -Inf(x);
  else {
         if (Sup(x) >= -Inf(x)) return Sup(x);
         else return -Inf(x);
       }
 
//  if (a > b)
//    return a;
//  else
//   return b;
}

/*----------------------------------------------------------------------------*/
REAL RDiam ( INTERVAL x )                                 // Relative diameter
{                                                         // of an interval
  if ( in(0.0,x) )                                        //------------------
    return Diam(x);
  else
    return( Diam(x) / AbsMin(x) );
}

/*----------------------------------------------------------------------------*/
REAL MaxRDiam ( INTERVAL_VECTOR& v )               // Maximum relative diameter
{                                                 //--------------------------
  REAL r;
  INT  i, l=1, u=Dimension(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RDiam(v(i)) > r) r = RDiam(v(i));
  return r;
}

//----------------------------------------------------------------------------
// Checks if the diameter of the interval 'x' is less or equal to 'n' ulps.
// Ulp is an abbreviation for units in the last place.
//----------------------------------------------------------------------------
INT UlpAcc ( INTERVAL x, INT n )
{
  INT  i;
  REAL Infimum;

  Infimum = Inf(x);
  for (i = 1; i <= n; i++) Infimum = Succ(Infimum);
  return (Infimum >= Sup(x));
}

//----------------------------------------------------------------------------
// Checks if the diameter of the interval vector 'v' is less or equal to 'n'
// ulps. Ulp is an abbreviation for: units in the last place of the mantissa.
//----------------------------------------------------------------------------
INT UlpAcc ( INTERVAL_VECTOR& v, INT n )
{
  INT  k, upper;

  for (k = 1, upper = Dimension(v); 
        (k < upper) && UlpAcc( v(k),n ); k++);
  return UlpAcc( v(k),n );
}

/*
interval Pi ( )                                    // Enclosure of constant pi
  { return(4.0*atan(_interval(1.0,1.0))); }        //-------------------------
*/

