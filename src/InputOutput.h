/*
 Purpose: Reading the problem and the command line, and writing the reports.

	ReadProblem(): Read an instance file and precompute the attraction
			terms of the existing facilities.
	GetArg(): Return the value following a command line option.
	ExistArg(): Determine whether a command line option is given.
	GetParams(): Process the command line and read the problem.
	ComputeResults(): Group the final boxes into the hulls of nearby boxes.
	OutPutStat(): Write the statistics of a run to a stream.

*******************************************************************************/

#ifndef __INPUTOUTPUT__
#define __INPUTOUTPUT__

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include "Matrix.h"
#include "Box.h"
#include "BLTree.h"
#include "Problems.h"

using namespace std;

PBLTREE ComputeResults( PBLTREE pbltv,INTERVAL_MATRIX & Results);

VOID OutPutStat (ostream & FileOut, BLTREEP_VECTOR & pbltv,
                 REAL Epsilon, int NFrac, int NDimDiv,
                 const VECTOR & BetterPoint, REAL BetterUbf,
                 double time, int NEvalFx, int NEvalFX, int NEvalG,
                 PCHAR NameProblem, bool MTest,
                 int NDelMPT, int NDelCTT, int NDelMT, int NDelFathGX);

PCHAR GetArg (PCHAR str_to_find, int narg, PPCHAR arg);
bool ExistArg (PCHAR str_to_find, int narg, PPCHAR arg);

void GetParams (int argc, PPCHAR argv, PCHAR NameProblem,  HTscalar_FctPtr & ff, 
                int & NDim, INTERVAL_VECTOR & XLim, VECTOR & Minima,
                REAL & Epsilon, int & NFrac, int & NDimDiv,
                bool & MTest, bool & EFatherGX, bool & UpdBau,
		int & WWidth, bool & FTest, bool & Draw);

#endif /*__INPUTOUTPUT__*/
