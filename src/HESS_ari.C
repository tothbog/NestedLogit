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
// File: HESS_ari (implementation)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to second order (i.e. Hessian matrix).
// Method: Overloading of operators and elementary functions for operations
//    of data types 'HessType' and 'HTvector'. To represent a Hessian matrix
//    by a lower triangular matrix, the type 'LowTriMatrix' is supplied.
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
#include <stdlib.h>         // For function 'exit()'
#include <iostream>       // I/O Handling   // .h was before
#include "Functions.h"      // Mathematical functions
#include "I_util.h"         // Interval matrix/vector utility functions
#include "HESS_ari.h"      

// The local variable 'HessOrder' is used to select the highest order of
// derivative which is computed. Its default value is 2, and normally
// the gradient and the Hessian matrix are computed.
//----------------------------------------------------------------------------
static INT HessOrder = 2;

//----------------------------------------------------------------------------
// Constructors, destructor, and resize function for the class 'LowTriMatrix'.
//----------------------------------------------------------------------------
VOID LowTriMatrix::init ( INT n )         // Initializes all components of the
{                                         // lower triangular matrix with 0.0
  if (n <= 0)                             //----------------------------------
    { dim = 0; rows = NULL; return; }
  dim = n;
  rows = new INTERVAL_VECTOR[n];
  for (INT i = 0; i < n; i++)                  // i-th row of the matrix takes
//    { Resize(rows[i],i+1); rows[i] = 0.0; }    // i components, i = 1,2,...,n

  { Resize(rows[i],i+1); 
    for (INT j = 1; j <= (i+1); j++) rows[i](j) = INTERVAL(0.0); }
     
}                                              //-----------------------------

LowTriMatrix::LowTriMatrix ( )                      // Default: no matrix rows
  { dim = 0; rows = NULL; }                         //------------------------


LowTriMatrix::LowTriMatrix ( INT n )             // Constructor: dimension n
  { init(n); }                              //-------------------------

LowTriMatrix::~LowTriMatrix ( )                                  // Destructor
  { delete [] rows; }                                            //-----------

VOID Resize ( LowTriMatrix& a, INT n )           // Resizes matrix to n rows
  { delete [] a.rows; a.init(n); }          //-------------------------

//----------------------------------------------------------------------------
// Assignment operators and copy constructor of the class 'LowTriMatrix'.
//----------------------------------------------------------------------------
LowTriMatrix& LowTriMatrix::operator= ( const LowTriMatrix& a )   // matrix = matrix
{                                                           //----------------
  if (this == &a) return *this;
  delete [] rows;
  dim = a.dim;
  if (dim == 0) { rows = NULL; return *this; }
  rows = new INTERVAL_VECTOR[dim];
  for (INT i = 0; i < dim; i++) rows[i] = a.rows[i];
  return *this;
}

LowTriMatrix& LowTriMatrix::operator= ( const REAL& r )       // matrix = real
{                                                             //--------------
  for (INT i = 0; i < dim; i++) {
    for (INT j = 1; j <= i+1; j++) 
      rows[i](j) = INTERVAL(r);
  }    
  return *this;
}

LowTriMatrix::LowTriMatrix ( LowTriMatrix& a )             // Copy constructor
{                                                          //-----------------
  rows = NULL;
  *this = a;   // use previously defined assignment operator
}

//----------------------------------------------------------------------------
// Component access for the type 'LowTriMatrix'
//----------------------------------------------------------------------------
INTERVAL_VECTOR& LowTriMatrix::operator() ( INT i )
{
  if ( (i < 1) || (i > dim) ) {
    cout << "Index out of range in "
         << "'INTERVAL_VECTOR& LowTriMatrix::operator() ( INT )'!" << endl;
    exit(-1);
  }
  return rows[i-1];      // Index range starts at 1
}

const INTERVAL_VECTOR& LowTriMatrix::operator() ( INT i ) const
{
  if ( (i < 1) || (i > dim) ) {
    cout << "Index out of range in "
         << "'ivector& LowTriMatrix::operator[] ( int )'!" << endl;
    exit(-1);
  }
  return rows[i-1];      // Index range starts at 1
}

//----------------------------------------------------------------------------
// Constructors, destructor, and resize function for the classes 'HessType'
// and 'HTvector'.
//----------------------------------------------------------------------------
HessType::HessType ( )                               // Default: no components
  { nmax = 0; }                                      //-----------------------

HessType::HessType ( INT n )                // Constructor: resizes internal
{                                             // arrays to dimension 'nmax'
  nmax = n;                            //------------------------------
  if (nmax <= 0) { nmax = 0; return; }
  if (HessOrder > 0)
     Resize(g,nmax);
  if (HessOrder > 1)
     Resize(h,nmax);
}

VOID Resize ( HessType& u, INT n )
{
  INT dim = n;
  u.nmax = dim;
  if (HessOrder > 0)
     Resize(u.g,dim);
  if (HessOrder > 1)
     Resize(u.h,dim);
}

HTvector::HTvector ( INT n )                   // Constructor: allocates 'n'
{                                                // 'HessType' components of
  nmax = n;                               // dimension 'n'
  if (nmax <= 0)                                 //---------------------------
    { nmax = 0; ht = NULL; return; }
  ht = new HessType[nmax];
  for (INT i = 0; i < nmax; i++) Resize(ht[i],nmax);
}

HTvector::~HTvector ( )                                          // Destructor
  { delete [] ht; }                                              //-----------

//----------------------------------------------------------------------------
// Assignment operators and copy constructors of the classes 'HessType' and
// 'HTvector'.
//----------------------------------------------------------------------------
HessType& HessType::operator= ( const HessType& u )           // HessType = HessType
{                                                       //--------------------
  if (this == &u) return *this;
  nmax = u.nmax;
  if (nmax > 0) { f = u.f;  g = u.g;  h = u.h; }
  return *this;
}

HessType::HessType ( const HessType& u )                         // Copy constructor
{                                                          //-----------------
  nmax = u.nmax;
  *this = u;   // use previously defined assignment operator
}

HessType& HessType::operator= ( const INTERVAL& x )          // HessType =  interval
{                                                      //---------------------
//  f = x;  g = 0.0;  h = 0.0;
//  return *this;
  
  f = x;
  if (HessOrder > 0)
    for (INT i = 1; i <= nmax; i++)  g(i) = INTERVAL(0.0);
  if (HessOrder > 1)
    h = 0.0;
  
  return *this;
}

HessType& HessType::operator= ( const REAL& r )            // HessType =  real
{                                                          //-----------------
  *this = INTERVAL(r);
  return *this;
}

HTvector& HTvector::operator= ( const HTvector& v )           // HTvector = HTvector
{                                                       //--------------------
  if (this == &v) return *this;
  delete [] ht;
  nmax = v.nmax;
  if (nmax == 0) { ht = NULL; return *this; }
  ht = new HessType[nmax];
  for (INT i = 0; i < nmax; i++) ht[i] = v.ht[i];
  return *this;
}

HTvector::HTvector ( const HTvector& v )                         // Copy constructor
{                                                          //-----------------
  ht = NULL;
  *this = v;   // use previously defined assignment operator
}

//----------------------------------------------------------------------------
// Component access for the type 'HTvector'.
//----------------------------------------------------------------------------
HessType& HTvector::operator() ( INT ind )
{
  INT n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'HessType& HTvector::operator() ( index )'!" << endl;
    exit(-1);
  }
  return ht[n-1];        // Index range starts at 1
}

const HessType& HTvector::operator() ( INT ind ) const
{
  int n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'HessType& HTvector::operator[] ( index )'!" << endl;
    exit(-1);
  }
  return ht[n-1];        // Index range starts at 1
}

//----------------------------------------------------------------------------
// Functions for range check of 'HessType' and 'HTvector' objects.
//----------------------------------------------------------------------------
VOID TestSize ( const HessType& a, const HessType& b, const char *fname )
{
  if (a.nmax != b.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

VOID TestSize ( const HTvector& v, const HTvector& w, const char *fname)
{
  if (v.nmax != w.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------------
// Transfer functions for variables
//----------------------------------------------------------------------------
HTvector HessVar ( const INTERVAL_VECTOR& x )                   // Generate variable
{                                                         //------------------
  INT       ubd = Dimension(x);
  HTvector  hht(ubd);

  for (INT i = 1; i <= ubd; i++) 
  {
    hht(i).f = x(i);
    if (HessOrder > 0)
      for (INT k = 1; k <= ubd; k++)
         hht(i).g(k) = (i == k) ? 1.0 : 0.0;
    if (HessOrder > 1)
    hht(i).h = 0.0;
  }
  return hht;
}

HTvector HessVar ( const VECTOR& v )                            // Generate variable
{                                                         //------------------
  INT     i, m = Dimension(v);
  INTERVAL_VECTOR u(m);

  for (i = 1; i <= m; i++) u(i) = v(i);
  return HessVar(u);
}

//----------------------------------------------------------------------------
// Access functions for the function value, the gradient, or the Hessian
//----------------------------------------------------------------------------
INTERVAL fValue ( const HessType& u )                          // Get function value
  { return u.f; }                                        //-------------------

INTERVAL_VECTOR gradValue ( const HessType& u )                        // Get gradient value
  { return u.g; }                                        //-------------------

INTERVAL_MATRIX hessValue ( const HessType& u )                         // Get Hessian value
{                                                         //------------------
  INTERVAL_MATRIX hm(u.nmax,u.nmax);

  for (INT i = 1; i <= u.nmax; i++) {
    hm(i,i) = u.h(i)(i);
    for (INT j=1; j <= i-1; j++) {
      hm(i,j) = u.h(i)(j);
      hm(j,i) = u.h(i)(j);
    }
  }
  return hm;
}

//----------------------------------------------------------------------------
// Access functions for vector-valued functions (function value or Jacobian)
//----------------------------------------------------------------------------
INTERVAL_VECTOR fValue ( const HTvector& u )                // Get function value of
{                                                     // n-dimensional vector-
  INTERVAL_VECTOR hv(u.nmax);                         // valued function
                                                      //----------------------
  for (INT i = 1; i <= u.nmax; i++) hv(i) = fValue(u(i));
  return hv;
}

INTERVAL_MATRIX JacValue ( const HTvector& u )              // Get Jacobian value of
{                                                     // n-dimensional scalar-
  INTERVAL_MATRIX hm(u.nmax,u.nmax);                  // valued function
                                                      //----------------------
  for (INT i = 1; i <= u.nmax; i++) 
//  hm[i] = gradValue(u[i]);
  SetRow( hm, i, gradValue(u(i)) );
  return hm;
}

//----------------------------------------------------------------------------
// Unary operators + and - for HessType operands
//----------------------------------------------------------------------------
HessType operator+ ( HessType& u )
  { return u; }

HessType operator- ( const HessType& u )
{
  HessType umin(u.nmax);

  umin.f = -u.f;
  if (HessOrder > 0)
    for (INT i = 1; i <= u.nmax; i++) {
      umin.g(i) = -u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) umin.h(i)(j) = -u.h(i)(j);
    }
  return umin;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for two HessType operands
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const HessType& v )
{
  HessType add(u.nmax);

  TestSize(u,v,"operator+ ( HessType&, HessType& )");
  add.f = u.f + v.f;
  if (HessOrder > 0)
    for (INT i = 1; i <= u.nmax; i++) {
      add.g(i) = u.g(i) + v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) add.h(i)(j) = u.h(i)(j) + v.h(i)(j);
    }
  return add;
}

HessType operator- ( const HessType& u, const HessType& v )
{
  HessType sub(u.nmax);

  TestSize(u,v,"operator- ( HessType&, HessType& )");
  sub.f = u.f - v.f;
  if (HessOrder > 0)
    for (INT i=1; i <= u.nmax; i++) {
      sub.g(i) = u.g(i) - v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) sub.h(i)(j) = u.h(i)(j) - v.h(i)(j);
    }
  return sub;
}

HessType operator* ( const HessType& u, const HessType& v )
{
  HessType mul(u.nmax);

  TestSize(u,v,"operator* ( HessType&, HessType& )");
  mul.f = u.f * v.f;
  if (HessOrder > 0)
    for (INT i = 1; i <= u.nmax; i++) {
      mul.g(i) = v.f*u.g(i) + u.f*v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          mul.h(i)(j) = v.f*u.h(i)(j) + u.g(i)*v.g(j) + v.g(i)*u.g(j)
                                                      + u.f*v.h(i)(j);
    }
  return mul;
}

HessType operator/ ( const HessType& u,const HessType& v )
{
  HessType div(u.nmax);

  TestSize(u,v,"operator/ ( HessType&, HessType& )");
  div.f = u.f / v.f;              // Can propagate 'division by zero' error
  if (HessOrder > 0)              //---------------------------------------
    for (INT i = 1; i <= u.nmax; i++) {
      div.g(i) = (u.g(i) - div.f*v.g(i)) / v.f;
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          div.h(i)(j) = (u.h(i)(j) - div.g(i)*v.g(j) - v.g(i)*div.g(j)
                                                     - div.f*v.h(i)(j)) / v.f;
    }
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for one interval and one HessType operand
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const INTERVAL& b )
{
  HessType add(u.nmax);

  add = u; add.f += b;
  return add;
}

HessType operator- ( const HessType& u, const INTERVAL& b )
{
  HessType sub(u.nmax);

  sub = u; sub.f -= b;
  return sub;
}

HessType operator* ( const HessType& u, const INTERVAL& b )
{
  HessType mul(u.nmax);

  mul.f = u.f * b;
  if (HessOrder > 0)
    for (INT i = 1; i <= u.nmax; i++) {
      mul.g(i) = b * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) mul.h(i)(j) = b * u.h(i)(j);
    }
  return mul;
}

HessType operator/ ( const HessType& u, const INTERVAL& b )
{
  HessType div(u.nmax);

  div.f = u.f / b;                // Can propagate 'division by zero' error
  if (HessOrder > 0)              //---------------------------------------
    for (INT i = 1; i <= u.nmax; i++) {
      div.g(i) = u.g(i) / b;
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) div.h(i)(j) = u.h(i)(j) / b;
    }
  return div;
}

HessType operator+ ( const INTERVAL& a, const HessType& v )
{
  HessType add(v.nmax);

  add = v; add.f = a + v.f;
  return add;
}

HessType operator- ( const INTERVAL& a, const HessType& v )
{
  HessType sub(v.nmax);

  sub.f = a - v.f;
  if (HessOrder > 0)
    for (INT i = 1; i <= v.nmax; i++) {
      sub.g(i) = -v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) sub.h(i)(j) = -v.h(i)(j);
    }
  return sub;
}

HessType operator* ( const INTERVAL& a, const HessType& v )
{
  HessType mul(v.nmax);

  mul.f = a * v.f;
  if (HessOrder > 0)
    for (INT i = 1; i <= v.nmax; i++) {
      mul.g(i) = a * v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++) mul.h(i)(j) = a * v.h(i)(j);
    }
  return mul;
}

HessType operator/ ( const INTERVAL& a, const HessType& v )
{
  HessType div(v.nmax);
  INTERVAL p, q;

  div.f = a / v.f;                // Can propagate 'division by zero' error
  if (HessOrder > 0) {            //---------------------------------------
    p = -div.f / v.f;
    q = (-2.0*p) / v.f;
    for (INT i = 1; i <= v.nmax; i++) {
      div.g(i) = p * v.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          div.h(i)(j) = p*v.h(i)(j) + q*v.g(i)*v.g(j);
    }
  }
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, * and / for one real and one HessType operand
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const REAL& b )
  { return u + INTERVAL(b); }

HessType operator- ( const HessType& u, const REAL& b )
  { return u - INTERVAL(b); }

HessType operator* ( const HessType& u, const REAL& b )
  { return u * INTERVAL(b); }

HessType operator/ ( const HessType& u, const REAL& b )
  { return u / INTERVAL(b); }     // Can propagate 'division by zero' error
                                  //---------------------------------------
HessType operator+ ( const REAL& a, const HessType& v )
  { return INTERVAL(a) + v; }

HessType operator- ( const REAL& a, const HessType& v )
  { return INTERVAL(a) - v; }

HessType operator* ( const REAL& a, const HessType& v )
  { return INTERVAL(a) * v; }

HessType operator/ ( const REAL& a, const HessType& v )
  { return INTERVAL(a) / v; }     // Can propagate 'division by zero' error
                                  //---------------------------------------

//----------------------------------------------------------------------------
// Elementary function for HessType arguments
//----------------------------------------------------------------------------
HessType Sqr ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1;

  res.f = Sqr(u.f);               // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 2.0 * u.f;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + 2.0*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Power ( const HessType& u, const INT k )
{
  HessType res(u.nmax);
  INTERVAL h1;

  if (k == 0)
    { res = 1.0; return res; }
  if (k == 1)
    return u;
  if (k == 2)
    return Sqr(u);

  res.f = Power(u.f,k);
  if (HessOrder > 0) {
    h1 = double(k) * Power(u.f,k-1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) +
                        double(k*(k-1))*Power(u.f,k-2)*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Power ( const HessType& u, const INTERVAL& k )
{
  HessType res(u.nmax);
  INTERVAL h1;
/*
  if (k == 0)
    { res = 1.0; return res; }
  if (k == 1)
    return u;
  if (k == 2)
    return Sqr(u);
*/
 
 res.f = Power(u.f,k);
  if (HessOrder > 0) {
    h1 = k * Power(u.f,k-1.0);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) +
                        k*(k-1.0)*Power(u.f,k-2.0)*u.g(i)*u.g(j);
    }
  }
  return res;
}


HessType Sqrt ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Sqrt(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = 0.5 / h0;
    h2 = -0.5*h1 / u.f;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Exp ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0;

  h0 = Exp(u.f);
  res.f = h0;
  if (HessOrder > 0)
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h0 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h0 * (u.h(i)(j) + u.g(i)*u.g(j));
    }
  return res;
}

HessType Log ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1,h2;

  res.f = Log(u.f);               // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/u.f;
    h2 = -Sqr(h1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Sin ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Sin(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = Cos(u.f);
    h2 = -h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Cos ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Cos(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = -Sin(u.f);
    h2 = -h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Tan ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Tan(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = Sqr(h0) + 1.0;                   // The subdistributive law implies
    h2 = 2.0*h0*h1;                       //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (INT i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g(i) = h1 * u.g(i);             //--------------------------------
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Cot ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Cot(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = -(Sqr(h0) + 1.0);                // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (INT i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g(i) = h1 * u.g(i);             //--------------------------------
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArcSin ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  res.f = ArcSin(u.f);            // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0-Sqr(u.f);
    h1 = 1.0/Sqrt(h0);
    h2 = u.f*h1/h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArcCos ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  res.f = ArcCos(u.f);            // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0-Sqr(u.f);
    h1 = -1.0/Sqrt(h0);
    h2 = u.f*h1/h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArcTan ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1, h2;

  res.f = ArcTan(u.f);
  if (HessOrder > 0) {
    h1 = 1.0/(1.0+Sqr(u.f));
    h2 = -2.0*u.f*Sqr(h1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArcCot ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1, h2;

  res.f = ArcCot(u.f);            // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = -1.0/(1.0+Sqr(u.f));
    h2 = 2.0*u.f*Sqr(h1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Sinh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Sinh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = Cosh(u.f);
    h2 = h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Cosh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Cosh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = Sinh(u.f);
    h2 = h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i)  = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Tanh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Tanh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = 1.0-Sqr(h0);                     // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (INT i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g(i)  = h1 * u.g(i);            //--------------------------------
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Coth ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  h0 = Coth(u.f);             // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = 1.0-Sqr(h0);                     // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (INT i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g(i) = h1 * u.g(i);             //--------------------------------
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArSinh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  res.f = ArSinh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0+Sqr(u.f);
    h1 = 1.0/Sqrt(h0);
    h2 = -u.f*h1/h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArCosh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h0, h1, h2;

  res.f = ArCosh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = Sqr(u.f)-1.0;
    h1 = 1.0/Sqrt(h0);
    h2 = -u.f*h1/h0;
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArTanh ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1, h2;

  res.f = ArTanh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/(1.0-Sqr(u.f));
    h2 = 2.0*u.f*Sqr(h1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType ArCoth ( const HessType& u )
{
  HessType res(u.nmax);
  INTERVAL h1,h2;

  res.f = ArCoth(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/(1.0-Sqr(u.f));
    h2 = 2.0*u.f*Sqr(h1);
    for (INT i = 1; i <= u.nmax; i++) {
      res.g(i) = h1 * u.g(i);
      if (HessOrder > 1)
        for (INT j = 1; j <= i; j++)
          res.h(i)(j) = h1*u.h(i)(j) + h2*u.g(i)*u.g(j);
    }
  }
  return res;
}

HessType Min ( const HessType& u, const HessType& v )
{
  HessType res(u.nmax);
  INTERVAL h1,h2;

  if (Sup(u.f) < Inf(v.f)) res.f = u.f;
  else if (Sup(v.f) < Inf(u.f)) res.f = v.f;
  else res.f = Hull( Inf(u.f) < Inf(v.f) ? Inf(u.f) : Inf(v.f) , Sup(u.f) < Sup(v.f) ? Sup(u.f) : Sup(v.f) );
  if (HessOrder > 0) {            //---------------------------
//	 res.g = u.g;
      if (Sup(u.f) < Inf(v.f)) res.g = u.g;
      else if (Sup(v.f) < Inf(u.f)) res.g = v.g;
      else res.g = Hull( u.g, v.g );
      if ( Inf(res.g) > Sup(res.g) ) 
      cerr << Inf(res.g) << endl << Sup(res.g) << endl;
}
      if (HessOrder > 1) {
	if (Sup(u.f) < Inf(v.f)) res.h = u.h;
	else if (Sup(v.f) < Inf(u.f)) res.h = v.h;
	else //res.g = Hull( u.h, v.h );
	  for (INT i = 1; i <= u.nmax; i++)
	    for (INT j = 1; j <= i; j++)
              res.h(i)(j) = Hull(u.h(i)(j), v.h(i)(j));

   } 
 
  return res;
}

// HessType Min ( const HessType& u, const INTERVAL & v )
// {
//   HessType res(u.nmax);
//   INTERVAL h1,h2;
// 
//   if (Sup(u.f) < Inf(v)) res.f = u.f;
//   else if (Sup(v) < Inf(u.f)) res.f = v;
//   else res.f = Hull( Inf(u.f) < Inf(v) ? Inf(u.f) : Inf(v) , Sup(u.f) < Sup(v) ? Sup(u.f) : Sup(v) );
//   if (HessOrder > 0) {            //---------------------------
// //	 res.g = u.g;
//       if (Sup(u.f) < Inf(v)) res.g = u.g;
//       else if (Sup(v) < Inf(u.f)) Clear(res.g) ;
//       else res.g = Hull( u.g, 0.0 );
//       if ( Inf(res.g) > Sup(res.g) ) 
//       cerr << Inf(res.g) << endl << Sup(res.g) << endl;
// }
//       if (HessOrder > 1) {
// 	if (Sup(u.f) < Inf(v)) res.h = u.h;
// 	else if (Sup(v) < Inf(u.f)) res.h=;
// 	else //res.g = Hull( u.h, v.h );
// 	  for (INT i = 1; i <= u.nmax; i++)
// 	    for (INT j = 1; j <= i; j++)
//               res.h(i)(j) = Hull(u.h(i)(j), 0.0);
// 
//    } 
//  
//   return res;
// }


HessType Max ( const HessType& u, const HessType& v )
{
  HessType res(u.nmax);
  INTERVAL h1,h2;

  if (Sup(u.f) < Inf(v.f)) res.f = v.f;
  else if (Sup(v.f) < Inf(u.f)) res.f = u.f;
  else res.f = Hull( Inf(u.f) > Inf(v.f) ? Inf(u.f) : Inf(v.f) , Sup(u.f) > Sup(v.f) ? Sup(u.f) : Sup(v.f) );
  if (HessOrder > 0) {            //---------------------------
//	 res.g = u.g;
      if (Sup(u.f) < Inf(v.f)) res.g = v.g;
      else if (Sup(v.f) < Inf(u.f)) res.g = u.g;
      else res.g = Hull( u.g, v.g );
      if ( Inf(res.g) > Sup(res.g) ) 
      cerr << Inf(res.g) << endl << Sup(res.g) << endl;
}
      if (HessOrder > 1) {
	if (Sup(u.f) < Inf(v.f)) res.h = v.h;
	else if (Sup(v.f) < Inf(u.f)) res.h = u.h;
	else //res.g = Hull( u.h, v.h );
	  for (INT i = 1; i <= u.nmax; i++)
	    for (INT j = 1; j <= i; j++)
              res.h(i)(j) = Hull(u.h(i)(j), v.h(i)(j));

   } 
 
  return res;
}


HessType ValueOrZero ( const HessType& u, const REAL& v )
{
  HessType res(u.nmax);
  INTERVAL h1,h2;

  if ( Sup(u.f) < v ) res.f = Hull(0.0);
  else if ( v <= Inf(u.f) ) res.f = u.f;
  else res.f = Hull( 0.0, Sup(u.f) );
  if (HessOrder > 0) {            //---------------------------
//	 res.g = u.g;
      if ( Sup(u.f) < v ) Clear(res.g); // = Hull(0.0);
      else if ( v < Inf(u.f) ) res.g = u.g;
      else for (int i=1; i<=u.nmax; i++) 
	res.g(i) = Hull( -1e255, 1e255 ); 
      if ( Inf(res.g) > Sup(res.g) ) 
      cerr << Inf(res.g) << endl << Sup(res.g) << endl;
}
      if (HessOrder > 1) {
	cerr << " Error: ValueOrZero with HessOrder > 1.\n";
	exit(0);
   } 
 
  return res;
}
//----------------------------------------------------------------------------
// Predefined routines for evaluation of HessType-functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
// Description: This function sets 'HessOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
VOID fEvalH ( HTscalar_FctPtr f, INTERVAL_VECTOR x, INTERVAL& fx )
{
  HessOrder = 0; 
  fx = fValue(f(HessVar(x))); 
  HessOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the gradient value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'gx': returns the gradient value 'grad f(x)'
// Description: This function sets 'HessOrder' to 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the gradient.
//----------------------------------------------------------------------------
VOID fgEvalH ( HTscalar_FctPtr f, INTERVAL_VECTOR x, INTERVAL& fx, 
                        	  INTERVAL_VECTOR& gx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  //  INT      n = Ub(x);
  INT	   n = Dimension(x);
  HessType fxH(n);

  HessOrder = 1;
  fxH       = f(HessVar(x));
  fx        = fValue(fxH);
  gx        = gradValue(fxH);
  HessOrder = 2;
}

//----------------------------------------------------------------------------
VOID fgEvalH1 ( HTscalar_FctPtr f, INTERVAL_VECTOR x, INTERVAL& fx, 
                        	  INTERVAL_VECTOR& gx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  //  INT      n = Ub(x);
  INT	   n = Dimension(x);
  HessType fxH(1);

  HessOrder = 1;
  fxH       = f(HessVar(x));
  fx        = fValue(fxH);
  gx        = gradValue(fxH);
  HessOrder = 2;
}
//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value, the gradient value, and the
//    Hessian matrix value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'gx': returns the gradient value 'grad f(x)'
//         'hx': returns the Hessian matrix value 'hess f(x)'
// Description: This function keeps 'HessOrder' = 2, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value, the
//    value of the gradient, and the value of the Hessian matrix.
//----------------------------------------------------------------------------
VOID fghEvalH ( HTscalar_FctPtr f, INTERVAL_VECTOR x, INTERVAL& fx, 
                                   INTERVAL_VECTOR& gx, INTERVAL_MATRIX& hx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  //  INT      n = Ub(x);
  INT      n = Dimension(x);
  HessType fxH(n);

  fxH = f(HessVar(x));
  fx  = fValue(fxH);
  gx  = gradValue(fxH);
  hx  = hessValue(fxH);
}

//----------------------------------------------------------------------------
VOID fghEvalH1 ( HTscalar_FctPtr f, INTERVAL_VECTOR x, INTERVAL& fx, 
                                   INTERVAL_VECTOR& gx, INTERVAL_MATRIX& hx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  //  INT      n = Ub(x);
  INT      n = Dimension(x);
  HessType fxH(1);

  fxH = f(HessVar(x));
  fx  = fValue(fxH);
  gx  = gradValue(fxH);
  hx  = hessValue(fxH);
}
//----------------------------------------------------------------------------
// Predefined routines for evaluation of HTvector-functions (with Jacobians)
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of type 'HTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the vector function value 'f(x)'
// Description: This function sets 'HessOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
VOID fEvalJ ( HTvector_FctPtr f, INTERVAL_VECTOR x, INTERVAL_VECTOR& fx )
{
  HessOrder = 0;
  fx        = fValue(f(HessVar(x)));
  HessOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the Jacobian matrix value.
// Parameters:
//    In : 'f' : function of type 'HTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'Jx': returns the Jacobian value 'Jf(x)'
// Description: This function sets 'HessOrder' to 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the Jacobian matrix.
//----------------------------------------------------------------------------
VOID fJEvalJ ( HTvector_FctPtr f, INTERVAL_VECTOR x, INTERVAL_VECTOR& fx, 
                                  INTERVAL_MATRIX& Jx )

{
  // HTvector fxGTv(Ub(x));     // not implemented in HP C++ V2.0!
  //  INT      n = Ub(x);
  INT      n = Dimension(x);
  HTvector fxGTv(n);

  HessOrder = 1;
  fxGTv     = f(HessVar(x));
  fx        = fValue(fxGTv);
  Jx        = JacValue(fxGTv);
  HessOrder = 2;
}

