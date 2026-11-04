#include "Problems.h"
#include "Bias1.h"
#include "BiasF.h"
#include "Functions.h"
#include "HESS_ari.h"
#include <iostream>

extern int n;                  // number of demand points
extern REAL b1, b2;            // parameters of the distance function
extern REAL a0, a1;            // parameters of the quality cost function G(alpha)
extern REAL c;                 // income per unit of captured demand (parameter F)
extern VECTOR w;               // buying power at the demand points
extern REAL ww;                // sum of the buying power
extern MATRIX fi;              // fi(i,2): parameter of the cost of the land G1(x) at demand point i
extern MATRIX p;               // locations of the demand points
extern VECTOR v;               // weight of the quality of the new facility at demand point i
extern VECTOR t2, r2;          // probabilistic rule: w(i)*attraction of the competitors,
                               // and attraction of all existing facilities
extern VECTOR aJ1, b;          // nested logit: sum of (attraction)^lambda_1 of the own chain,
                               // and the aggregated attraction of the competing chains
extern REAL la1;               // lambda_1 of the nested logit model
extern REAL Beta;              // distance decay: g(d) = exp(Beta*d)

/*----------------------------------------------------------------------------*/
/* Nested logit choice rule.                                                   */
/*                                                                             */
/*   F(x) = G1(x1,x2) + G2(alpha) - c * M(x)                                   */
/*                                                                             */
/* where M(x) is the demand captured by the locating chain:                    */
/*   M(x) = sum_i w(i) * (1 - b(i) / (A(i)^(1/lambda_1) + b(i))),               */
/*   A(i) = (v(i)*alpha/g(d_i(x)))^lambda_1 + aJ1(i),                           */
/* G1(x) = sum_i w(i)/(d_i(x)^2 + fi(i,2)) is the cost of the land and          */
/* G2(alpha) = exp(alpha/a0 + a1) - exp(a1) is the cost of quality.             */
/*----------------------------------------------------------------------------*/
HessType fFacLoc(const HTvector &x) {
  HessType psisum(3), Msum(3), hh(3);
  HessType dd(3), gdd(3);
  INTERVAL La1=Hull(la1), invLa1 = 1.0 / La1;
  int i;

  psisum = 0.0;
  Msum = 0.0;
  for (i = 1; i <= n; i++) {
    dd = b1 * Sqr(x(1) - p(i, 1)) + b2 * Sqr(x(2) - p(i, 2)); // squared distance

    psisum = psisum + w(i) / (dd + fi(i, 2));                 // cost of the land

    gdd = Exp(Beta * Sqrt(dd + 1e-255));                      // g(d) = e^(Beta*d)

    // demand of point i lost to the competing chains
    Msum = Msum +
           w(i) * b(i) /
               (Power(Power(v(i) * x(3) / gdd, La1) + aJ1(i), invLa1) + b(i));
  }

  hh = psisum + (Exp(x(3) / a0 + a1) - exp(a1)) - c * ww + c * Msum;
  return hh;
}

/*----------------------------------------------------------------------------*/
/* Probabilistic choice rule: the new facility competes with all existing     */
/* facilities in a single step, chains are not taken into account.            */
/*----------------------------------------------------------------------------*/
HessType fFacLocViejo(const HTvector &x) {
  HessType psisum(3), Msum(3), hh(3);
  HessType dd(3), gdd(3);
  int i;

  psisum = 0.0;
  Msum = 0.0;
  for (i = 1; i <= n; i++) {
    dd = b1 * Sqr(x(1) - p(i, 1)) + b2 * Sqr(x(2) - p(i, 2));

    psisum = psisum + w(i) / (dd + fi(i, 2));

    gdd = Exp(Sqrt(dd + 1e-255));

    Msum = Msum + t2(i) / (v(i) * x(3) / gdd + r2(i));
  }

  hh = psisum + (Exp(x(3) / a0 + a1) - exp(a1)) - c * ww + c * Msum;
  return hh;
}

/*----------------------------------------------------------------------------*/
/* fFacLoc, printing the parts of the objective.                               */
/*----------------------------------------------------------------------------*/
HessType fFacLocD(const HTvector &x)
{
  HessType psisum(3), Msum(3), hh(3);
  HessType dd(3), gdd(3);
  INTERVAL La1=Hull(la1), invLa1 = 1.0 / La1;
  int i;

  psisum = 0.0;
  Msum = 0.0;
  for (i = 1; i <= n; i++) {
    dd = b1 * Sqr(x(1) - p(i, 1)) + b2 * Sqr(x(2) - p(i, 2));

    psisum = psisum + w(i) / (dd + fi(i, 2));

    gdd = Exp(Beta * Sqrt(dd + 1e-255));

    Msum =
        Msum +
        w(i) * b(i) /
            (Power(Power(v(i) * x(3) / gdd, La1) + aJ1(i), invLa1) + b(i));
  }

  cerr << "c*M = [ " << Inf(fValue(c * ww - c * Msum)) << " , "
       << Sup(fValue(c * ww - c * Msum)) << " ]" << endl;
  cerr << "G1(X) = [ " << Inf(fValue(psisum)) << " , " << Sup(fValue(psisum))
       << " ]" << endl;
  cerr << "G2(alpha) = [ " << Inf(fValue(Exp(x(3) / a0 + a1) - exp(a1)))
       << " , " << Sup(fValue(Exp(x(3) / a0 + a1) - exp(a1))) << " ]" << endl;

  hh = psisum + (Exp(x(3) / a0 + a1) - exp(a1)) - c * ww + c * Msum;

  cerr << "F_NestedLogit = [ " << Inf(fValue(hh)) << " , " << Sup(fValue(hh))
       << " ]" << endl;
  return hh;
}

/*----------------------------------------------------------------------------*/
/* fFacLocViejo, printing the parts of the objective.                          */
/*----------------------------------------------------------------------------*/
HessType fFacLocViejoD(const HTvector &x) {
  HessType psisum(3), Msum(3), hh(3);
  HessType dd(3);
  int i;

  psisum = 0.0;
  Msum = 0.0;
  for (i = 1; i <= n; i++) {
    dd = b1 * Sqr(x(1) - p(i, 1)) + b2 * Sqr(x(2) - p(i, 2));

    psisum = psisum + w(i) / (dd + fi(i, 2));

    dd = Exp(Sqrt(dd + 1e-255));

    Msum = Msum + t2(i) / (v(i) * x(3) / dd + r2(i));
  }
  cerr << "c*M_Prob = [ " << Inf(fValue(c * ww - c * Msum)) << " , "
       << Sup(fValue(c * ww - c * Msum)) << " ]" << endl;
  hh = psisum + (Exp(x(3) / a0 + a1) - exp(a1)) - c * ww + c * Msum;
  cerr << "F_Prob = [ " << Inf(fValue(hh)) << " , " << Sup(fValue(hh)) << " ]"
       << endl;

  return hh;
}
