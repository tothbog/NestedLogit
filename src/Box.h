/*
 Box.h - the BOX container of the branch-and-bound algorithm.

 A BOX stores
   X          : the subregion of the search space (x1, x2, alpha),
   FX         : an inclusion of the objective function over X,
   M          : flag set when the box was found monotone,
   feas       : true if the whole box is known to be feasible,
   fe         : fe(i) is true if the box is known to satisfy the i-th
                (demand point) feasibility constraint,
   NBox       : identifier of the box (used in the graphical output),
   FatherNBox : identifier of the box it was obtained from.

 EvalpF() computes the pf value (BetterUbf - lb F(X)) / w(F(X)) of a box.
 The Draw()/UnDraw() members write the box for the graphical output
 (the -tcl option and the *_res.tclin file).
*/

#ifndef __BOX__
#define __BOX__

#include "IntervalVector.h"
#include "IntervalVectorPlusz.h"
#include "IntegerVector.h"
#include "Matrix.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include <stdio.h>

typedef char * PCHAR; 
typedef PCHAR * PPCHAR; 

class BOX;

typedef BOX *	PBOX; 
typedef PBOX *	PPBOX;

class BOX 
{
 public:
   INTERVAL_VECTOR X;
   INTERVAL FX;
   bool M;
   bool feas;
   INT NBox;
   INT FatherNBox;
   INTEGER_VECTOR fe;

 public:
   BOX (INT nn);
   BOX (const INTERVAL_VECTOR & IV);
   BOX (const BOX & B);
   BOX (bool son, const BOX & B);   // copy; if son, it gets a new NBox
   ~BOX ();

   friend int Dimension  (const BOX & );
   friend REAL EvalpF	 (const BOX & B, int Index, REAL BetterUbf);
   void Draw      (const INTERVAL_VECTOR & XLim, INT WWidth, PCHAR color);
   void Draw      (const INTERVAL_VECTOR & XLim, INT WWidth, PCHAR color, ostream & out);
   void UnDraw    ();
   friend ostream & operator << (ostream & os, const BOX & B);
   friend ostream & operator << (ostream & os, PBOX & pb );   
};

#endif /*__BOX__*/
