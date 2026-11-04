/*
 Purpose: Objective functions of the facility location and design problem,
	evaluated with the automatic differentiation of the C++ Toolbox for
	Verified Computing (HessType arithmetic, adapted to PROFIL/BIAS).

 The variables are x = (x1, x2, alpha): location and quality of the new
 facility. All functions return the negative of the profit of the locating
 chain, so the problem is solved as a minimisation problem.

 fFacLoc():       nested logit choice rule (the model of the article).
 fFacLocViejo():  probabilistic (Huff-like) choice rule, for comparison (-vi).
 fFacLocD(), fFacLocViejoD(): the same functions, also printing the terms of
                  the objective (used for the final report and for -ft).

*******************************************************************************/

#ifndef __PROBLEMS__
#define __PROBLEMS__

#include "Box.h"
#include "HESS_ari.h"

/*----------------------------------------------------------------------------*/
HessType fFacLoc        ( const HTvector& x );
HessType fFacLocViejo   ( const HTvector& x );
HessType fFacLocD       ( const HTvector& x );
HessType fFacLocViejoD  ( const HTvector& x );

#endif /*__PROBLEMS__*/
