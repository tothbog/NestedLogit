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
// File: XI_ari (implementation)
// Purpose: Definition of an extended interval arithmetic which allows the
//    division by an interval containing zero.
// Method: Overloading of operators for arithmetic and lattice operations
//    of data type 'XINTERVAL'.
// Global type:
//    KindType         : component type of extended intervals
// Global function:
//    EmptyIntval()    : empty set as irregular interval
// Class XINTERVAL:
//    operators %, -, &: operators of extended interval arithmetic
//    operator =       : assignment operator
//----------------------------------------------------------------------------
#include "I_util.h"     // Interval utility functions
#include "XI_ari.h"

//----------------------------------------------------------------------------
// An extended interval 'x', represented by the type 'xinterval', is defined
// according to the following rules (a <= b):
//
// x = [a, b]              :  x.kind = Finite,     x.inf = a, x.sup = b
// x = [a, +oo]            :  x.kind = PlusInfty,  x.inf = a, x.sup undefined
// x = [-oo, a]            :  x.kind = MinusInfty, x.sup = a, x.inf undefined
// x = [-oo, a] v [b, +oo] :  x.kind = Double,     x.inf = b, x.sup = a
// x = [-oo, +oo]          :  x.kind = Double,     x.inf = a, x.sup = a
//
// In this definition, 'v' stands for the set union and 'oo' for infinity.
//----------------------------------------------------------------------------
INTERVAL EmptyIntval ( )                         // Irregular (empty) interval
{                                                //---------------------------
  INTERVAL x;
//x = _UNCHECKED_INTERVAL(999999999.0,-999999999.0);
  x = INTERVAL(999999999.0,-999999999.0);
  return x;
}

//----------------------------------------------------------------------------
// Constructors and assignment operator
//----------------------------------------------------------------------------
XINTERVAL::XINTERVAL ( )
{
  kind = Finite;
  inf  = 0.0;
  sup  = 0.0;
}

XINTERVAL::XINTERVAL ( const KindType& k, const REAL& i, const REAL& s )
{
  kind = k;
  inf  = i;
  sup  = s;
}

XINTERVAL::XINTERVAL ( const XINTERVAL& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
}

XINTERVAL& XINTERVAL::operator= ( const XINTERVAL& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
  return *this;
}

//----------------------------------------------------------------------------
// Subtraction of an extended interval 'B' from a real value 'a'.
//----------------------------------------------------------------------------
XINTERVAL operator- ( const REAL& a, const XINTERVAL& B )
{
  XINTERVAL D;

  switch (B.kind) {
//    case Finite     : D.kind = Finite;                    // D = [D.inf,D.sup]	 
//                      D.inf  = subd(a,B.sup);             //------------------
//                      D.sup  = subu(a,B.inf);
//                      break;
    case Finite     : D.kind = Finite;                    // D = [D.inf,D.sup]	 
                      BiasRoundDown();
                      D.inf  = a - B.sup;             
                      BiasRoundUp();
                      D.sup  = a - B.inf;
                      BiasRoundNear();
                      break;
   
    case PlusInfty  : D.kind = MinusInfty;                    // D = [inf,+oo]
                      BiasRoundUp();
                      D.sup  = a - B.inf;                 
                      BiasRoundNear();
                      break;
    case MinusInfty : D.kind = PlusInfty;                     // D = [-oo,sup]
                      BiasRoundDown();
                      D.inf  = a - B.sup;                
                      BiasRoundNear();
                      break;
    case Double    :  D.kind = Double;        // D = [-oo,D.sup] v [D.inf,+oo]
                      BiasRoundDown();
                      D.inf  = a - B.sup; 
                      BiasRoundUp();
                      D.sup  = a - B.inf;
                      BiasRoundNear();
                      if (D.inf < D.sup) D.inf = D.sup;
                      break;
  } // switch
  return D;
}

//----------------------------------------------------------------------------
// Extended interval division 'A / B' where 0 in 'B' is allowed.
//----------------------------------------------------------------------------
XINTERVAL operator% ( const INTERVAL& A, const INTERVAL& B )
{
  INTERVAL  c;
  XINTERVAL Q;

  if ( in(0.0,B) ) {
    if ( ( (Inf(A) < 0.0) && (0.0 < Sup(A) ) ) || (A == 0.0) || (B == 0.0) ) {
      Q.kind = Double;                    // Q = [-oo,+oo] = [-oo,0] v [0,+oo]
      Q.sup  = 0.0;                       //----------------------------------
      Q.inf  = 0.0;
    }
    else if ( (Sup(A) <= 0.0) && (Sup(B) == 0.0) ) {        // Q = [Q.inf,+oo]
      Q.kind = PlusInfty;                                   //----------------
      BiasRoundDown();
      Q.inf  = Sup(A) / Inf(B);
      BiasRoundNear();
    }
    else if ( (Sup(A) <= 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      BiasRoundUp();
      Q.sup  = Sup(A) / Sup(B);        
      BiasRoundDown();
      Q.inf  = Sup(A) / Inf(B);
      BiasRoundNear();
    }
    else if ( (Sup(A) <= 0.0) && (Inf(B) == 0.0) ) {        // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      BiasRoundUp();
      Q.sup  = Sup(A) / Sup(B);
      BiasRoundNear();
    }
    else if ( (Inf(A) >= 0.0) && (Sup(B) == 0.0) ) {        // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      BiasRoundUp();
      Q.sup  = Inf(A) / Inf(B);
      BiasRoundNear();
    }
    else if ( (Inf(A) >= 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      BiasRoundUp();
      Q.sup  = Inf(A) / Inf(B);
      BiasRoundDown();
      Q.inf  = Inf(A) / Sup(B);
      BiasRoundNear();
    }
    else { // if ( (Inf(A) >= 0.0) && (Inf(B) == 0.0) )
      Q.kind = PlusInfty;                                   // Q = [Q.inf,+oo]
      BiasRoundDown();
      Q.inf  = Inf(A) / Sup(B);
      BiasRoundNear();
    }
  } // in(0.0,B)
  else {  // !in(0.0,B)
    c = A / B;                                            // Q = [C.inf,C.sup]
    Q.kind = Finite;                                      //------------------
    Q.inf  = Inf(c);
    Q.sup  = Sup(c);
  }

  return Q;
} // operator%

//----------------------------------------------------------------------------
// Intersection of an interval 'X' and an extended interval 'Y'. The result
// is given as a pair (vector) of intervals, where one or both of them can
// be empty intervals.
//----------------------------------------------------------------------------
INTERVAL_VECTOR operator& ( const INTERVAL& X, const XINTERVAL& Y )
{
  INTERVAL         H, ISECT;
  INTERVAL_VECTOR  IS(2);

  IS(1) = EmptyIntval();
  IS(2) = EmptyIntval();

  switch (Y.kind) {
    case Finite     : // [X.inf,X.sup] & [Y.inf,Y.sup]
                      //------------------------------
                      H = INTERVAL(Y.inf,Y.sup);
                   //   if ( !Disjoint(X,H) ) 
                   //   IS(1) = X & H;
                   //   break;
                      if (Intersection(ISECT,X,H)==1)   // the same as
                       IS(1) = ISECT;                   // !Disjoint(X,H)
                      break;
    case PlusInfty  : // [X.inf,X.sup] & [Y.inf,+oo]
                      //----------------------------
                      if (Sup(X) >= Y.inf)
                        if (Inf(X) > Y.inf)
                          IS(1) = X;
                        else
                          IS(1) = INTERVAL(Y.inf,Sup(X));
                      break;
    case MinusInfty : // [X.inf,X.sup] & [-oo,Y.sup]
                      //----------------------------
                      if (Y.sup >= Inf(X))
                        if (Sup(X)<Y.sup)
                          IS(1) = X;
                        else
                          IS(1) = INTERVAL(Inf(X),Y.sup);
                      break;
    case Double     : if ( (Inf(X) <= Y.sup) && (Y.inf <= Sup(X)) ) {
                        IS(1) = INTERVAL(Inf(X),Y.sup);    // X & [-oo,Y.sup]
                        IS(2) = INTERVAL(Y.inf,Sup(X));    // X & [Y.inf,+oo]
                      }
                      else if (Y.inf <= Sup(X)) // [X.inf,X.sup] & [Y.inf,+oo]
                        if (Inf(X) >= Y.inf)    //----------------------------
                          IS(1) = X;
                        else
                          IS(1) = INTERVAL(Y.inf,Sup(X));
                      else if (Inf(X) <= Y.sup) // [X.inf,X.sup] & [-oo,Y.sup]
                        if (Sup(X) <= Y.sup)    //----------------------------
                          IS(1) = X;
                        else
                          IS(1) = INTERVAL(Inf(X),Y.sup);
                      break;
  } // switch

  return IS;
} // operator&

