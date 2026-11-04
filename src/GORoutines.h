/*
 Purpose: Common routines of the interval branch-and-bound algorithm.

 InitBox(), InitBoxGrad(): Set the initial box to the search region and
        evaluate the objective (and its gradient) on it.
 BaumannPoint(): Compute the Baumann point of a box from its gradient
        inclusion (used with the -ub option and on the final boxes).
 UpdateBetterUbf(): If the objective value at the point pc is better than the
        current upper bound of the global minimum, update it.
 CutOffTest(): Remove from a tree all boxes which can not contain a global
        minimiser. Called after the upper bound is improved.
 CutOffTheActualBox(): Reject the current box if its lower bound is above the
        upper bound of the global minimum.
 MonoTestFeas(): Monotonicity test. If the objective is monotone in some
        variable over the box and the best facet is feasible, the box is
        deleted. If the box lies on the border of the search region, it is
        reduced to the border facet instead.
 Feasible(): Check whether the box is certainly feasible (1), certainly
        infeasible (-1) or undetermined (0) with respect to the forbidden
        circles around the demand points.

*******************************************************************************/

#ifndef __GOROUTINES__
#define __GOROUTINES__

#include "Vector.h"
#include "Box.h"
#include "BLTree.h"
#include "Problems.h"
#include "XI_ari.h"
#include "float.h"
#include "math.h"

/*----------------------------------------------------------------------------*/
void InitBox ( PBOX & pb, const INTERVAL_VECTOR & XLim, HTscalar_FctPtr ff,
               int & NEvalFX );
void InitBoxGrad ( PBOX & pb, INTERVAL_VECTOR & GX, const INTERVAL_VECTOR & XLim,
		HTscalar_FctPtr ff, int & NEvalG );
void BaumannPoint(const PBOX & pb,INTERVAL_VECTOR & GX,  PBOX & pc);
bool UpdateBetterUbf (const PBOX & pb, REAL & BetterUbf, VECTOR & BetterPoint, 
                      HTscalar_FctPtr ff, int & NEvalFx, PBOX & pc);
void CutOffTest (PBLTREE pblt, REAL BetterUbf, int & NDelCTT, 
                 const INTERVAL_VECTOR & XLim, INT WWidth, 
                 bool Draw, PCHAR color );
bool CutOffTheActualBox ( PBOX & pb, REAL BetterUbf, int & NDelMPT,
                    const INTERVAL_VECTOR & XLim, INT WWidth, bool Draw);
bool MonoTestFeas(PBOX & pb, INTERVAL_VECTOR & GX, int & NDelMT, const INTERVAL_VECTOR & XLim, 
			INT WWidth, bool Draw);
int Feasible (PBOX & pb);

#endif /*__GOROUTINES__*/
