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
// File: HESS_ari (header)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to second order (i.e. Hessian matrix).
// Types:
//    HTscalar_FctPtr         : pointer for a scalar valued function
//    HTvector_FctPtr         : pointer for a vector valued function
// Class LowTriMatrix:
//    LowTriMatrix()          : constructors
//    ~LowTriMatrix()         : destructor
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for triangular matrices
//                              and real values
//    operator ()             : component access
// Class HessType:
//    HessType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              HessType, interval, and real
//    Dim()                   : returns the dimension of the Hessian matrix
//    HessVar()               : to define HessType variables
//    fValue(), gradValue(),
//    hessValue()             : to get function and derivative values
//    operators +, -, *, /    : operators of diff. arithmetic
//    Sqr(), Sqrt(), Power(),
//    Exp(), Sin(), Cos(), ...: elementary functions of diff. arithmetic
//    fEvalH()                : to compute function value only
//    fgEvalH()               : to compute function and first derivative
//                              value (gradient)
//    fghEvalH()              : to compute function value, gradient, and
//                              Hessian matrix value
// Class HTvector:
//    HTvector()              : constructors
//    ~HTvector()             : destructor
//    operator =              : assignment operator
//    operator ()             : component access
//    Dim()                   : to get the actual dimension
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian)
//----------------------------------------------------------------------------
#ifndef __HESS_ARI__
#define __HESS_ARI__

#include "Interval.h"         // Interval arithmetic
#include "IntervalMatrix.h"   // Interval matrix/vector arithmetic
#include "I_util.h"           // Interval matrix/vector utility functions
#include "Global.h"           // Class 'index'

class LowTriMatrix {
    INT             dim;
    INTERVAL_VECTOR *rows;

    VOID init ( INT );
    friend VOID Resize ( LowTriMatrix&, INT );

  public:
    LowTriMatrix ( );
    explicit LowTriMatrix ( INT );
    LowTriMatrix ( LowTriMatrix& );
    ~LowTriMatrix ( );

    LowTriMatrix&    operator= ( const LowTriMatrix& );
    LowTriMatrix&    operator= ( const REAL& );
    INTERVAL_VECTOR& operator() ( INT );
    const INTERVAL_VECTOR& operator() ( INT ) const;

};

class HessType;
class HTvector;

typedef HessType (*HTscalar_FctPtr)(const HTvector&);
typedef HTvector (*HTvector_FctPtr)(const HTvector&);

class HessType {
  
    INT              nmax;
    INTERVAL         f;
    INTERVAL_VECTOR  g;
    LowTriMatrix     h;

    friend VOID Resize ( HessType&, INT );
    friend VOID TestSize ( const HessType&, const HessType&, const char* );

  public:
    HessType ( );
    HessType ( INT );
    HessType ( const HessType& );

    HessType& operator= ( const HessType& );
    HessType& operator= ( const INTERVAL& );
    HessType& operator= ( const REAL& );

    INT  Dim ( ) { return nmax; };

    friend HTvector HessVar ( const INTERVAL_VECTOR& );
    friend HTvector HessVar ( const VECTOR& );

    friend INTERVAL fValue    ( const HessType& );
    friend INTERVAL_VECTOR  gradValue ( const HessType& );
    friend INTERVAL_MATRIX  hessValue ( const HessType& );

    friend HessType operator+ ( HessType&);
    friend HessType operator- ( const HessType&);
    friend HessType operator+ ( const HessType&, const HessType&);
    friend HessType operator- ( const HessType&, const HessType&);
    friend HessType operator* ( const HessType&, const HessType&);
    friend HessType operator/ ( const HessType&, const HessType&);
    friend HessType operator+ ( const HessType&, const INTERVAL&);
    friend HessType operator- ( const HessType&, const INTERVAL&);
    friend HessType operator* ( const HessType&, const INTERVAL&);
    friend HessType operator/ ( const HessType&, const INTERVAL&);
    friend HessType operator+ ( const INTERVAL&, const HessType&);
    friend HessType operator- ( const INTERVAL&, const HessType&);
    friend HessType operator* ( const INTERVAL&, const HessType&);
    friend HessType operator/ ( const INTERVAL&, const HessType&);
    friend HessType operator+ ( const HessType&, const REAL&);
    friend HessType operator- ( const HessType&, const REAL&);
    friend HessType operator* ( const HessType&, const REAL&);
    friend HessType operator/ ( const HessType&, const REAL&);
    friend HessType operator+ ( const REAL&, const HessType&);
    friend HessType operator- ( const REAL&, const HessType&);
    friend HessType operator* ( const REAL&, const HessType&);
    friend HessType operator/ ( const REAL&, const HessType&);

    friend HessType Sqr    ( const HessType& );
    friend HessType Power  ( const HessType&, const INT );
    friend HessType Power  ( const HessType&, const INTERVAL& );
    friend HessType Sqrt   ( const HessType& );
    friend HessType Exp    ( const HessType& );
    friend HessType Log    ( const HessType& );
    friend HessType Sin    ( const HessType& );
    friend HessType Cos    ( const HessType& );
    friend HessType Tan    ( const HessType& );
    friend HessType Cot    ( const HessType& );
    friend HessType ArcSin ( const HessType& );
    friend HessType ArcCos ( const HessType& );
    friend HessType ArcTan ( const HessType& );
    friend HessType ArcCot ( const HessType& );
    friend HessType Sinh   ( const HessType& );
    friend HessType Cosh   ( const HessType& );
    friend HessType Tanh   ( const HessType& );
    friend HessType Coth   ( const HessType& );
    friend HessType ArSinh ( const HessType& );
    friend HessType ArCosh ( const HessType& );
    friend HessType ArTanh ( const HessType& );
    friend HessType ArCoth ( const HessType& );
    friend HessType Min    ( const HessType&, const HessType& );
//    friend HessType Min ( const HessType& u, const INTERVAL & v );
    friend HessType Max    ( const HessType&, const HessType& );
    friend HessType ValueOrZero ( const HessType& u, const REAL& v );
    
    friend VOID fEvalH   ( HTscalar_FctPtr, INTERVAL_VECTOR, INTERVAL& );
    friend VOID fgEvalH  ( HTscalar_FctPtr, INTERVAL_VECTOR, INTERVAL&, 
                           INTERVAL_VECTOR& );
    friend VOID fgEvalH1  ( HTscalar_FctPtr, INTERVAL_VECTOR, INTERVAL&, 
                           INTERVAL_VECTOR& );
    friend VOID fghEvalH ( HTscalar_FctPtr, INTERVAL_VECTOR, INTERVAL&, 
                           INTERVAL_VECTOR&, INTERVAL_MATRIX& );
    friend VOID fghEvalH1 ( HTscalar_FctPtr, INTERVAL_VECTOR, INTERVAL&, 
                           INTERVAL_VECTOR&, INTERVAL_MATRIX& );
};

class HTvector {
    INT       nmax;
    HessType *ht;

    friend VOID TestSize ( const HTvector&, const HTvector&, const char* );

  public:
    HTvector  ( INT );
    HTvector  ( const HTvector& );
    ~HTvector ( );

    HTvector& operator=  ( const HTvector& );
    HessType& operator() ( INT );
    const HessType& operator() ( INT ) const;
    
    INT Dim ( ) { return nmax; };

    friend INTERVAL_VECTOR  fValue   ( const HTvector& );
    friend INTERVAL_MATRIX  JacValue ( const HTvector& );

    friend VOID fEvalJ  ( HTvector_FctPtr, INTERVAL_VECTOR, 
                          INTERVAL_VECTOR& );                               
    friend VOID fJEvalJ ( HTvector_FctPtr, INTERVAL_VECTOR, 
                          INTERVAL_VECTOR&, INTERVAL_MATRIX& );
};
#endif
