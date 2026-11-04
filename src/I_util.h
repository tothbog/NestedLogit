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
//----------------------------------------------------------------------------
// File: I_util (header)
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
#ifndef __I_UTIL__
#define __I_UTIL__

#include "Interval.h"
#include "Vector.h"
#include "IntervalVector.h"

bool             in	      ( REAL, const INTERVAL & );
bool             in	      ( const INTERVAL &, const INTERVAL & );
bool             in       ( const INTERVAL_VECTOR &, const INTERVAL_VECTOR &);
bool             in       ( INT, const INTERVAL_VECTOR & );
bool			 in 	  ( const VECTOR &, const INTERVAL_VECTOR &);
INTERVAL         Blow     ( INTERVAL, REAL);
INTERVAL_VECTOR  Blow     ( INTERVAL_VECTOR&, REAL );
INT	             Disjoint ( INTERVAL, INTERVAL);
INT              Disjoint ( INTERVAL_VECTOR&, INTERVAL_VECTOR& );
REAL	         AbsMin   ( INTERVAL );
REAL	  	     AbsMax   ( INTERVAL );
REAL	 	 	 RDiam    ( INTERVAL );
REAL     	 	 MaxRDiam ( INTERVAL_VECTOR& );
INT 	  	 	 UlpAcc   ( INTERVAL, INT );
INT      	 	 UlpAcc   ( INTERVAL_VECTOR&, INT );

#endif
