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
// File: XI_ari (header)
// Purpose: Definition of an extended interval arithmetic which allows the
//    division by an interval containing zero.
// Global type:
//    KindType         : component type of extended intervals
// Global function:
//    EmptyIntval()    : empty set as irregular interval
// Class XINTERVAL:
//    operators %, -, &: operators of extended interval arithmetic
//    operator =       : assignment operator
//----------------------------------------------------------------------------
#ifndef __XI_ARI__
#define __XI_ARI__

#include "Interval.h"         // Interval arithmetic
#include "IntervalVector.h"   // Interval vector arithmetic

typedef enum { Finite, PlusInfty, MinusInfty, Double } KindType;

INTERVAL EmptyIntval ( );   // Irregular (empty) INTERVAL

class XINTERVAL {                  // Extended INTERVALs according
  private:                         // to the above definition
    KindType  kind;                //-----------------------------
    REAL      inf, sup;

  public:
    XINTERVAL ( );
    XINTERVAL ( const KindType&, const REAL&, const REAL& );
    XINTERVAL ( const XINTERVAL& );

    XINTERVAL& operator= ( const XINTERVAL& );

    friend XINTERVAL         operator- ( const REAL&, const XINTERVAL& );
    friend XINTERVAL         operator% ( const INTERVAL&, const INTERVAL& );
    friend INTERVAL_VECTOR   operator& ( const INTERVAL&, const XINTERVAL& );
};
#endif
