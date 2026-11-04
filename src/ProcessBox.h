/*
 Purpose: Processing of a box: rejection tests, division, and storing it in
	the working or in the final data structure.

 CheckBox(): Apply the rejection tests to a new box and, if it is not
	rejected, store it in the working tree or, if it satisfies the
	termination criterion, in the final tree.
 DivideBox(): Divide the current box into subboxes and check each of them
	with CheckBox(). NDimDiv is the number of (widest) dimensions to divide,
	NFrac is the number of parts in each divided dimension.

*******************************************************************************/

#ifndef __PROCESSBOX__
#define __PROCESSBOX__

#include "float.h"
#include "Box.h"
#include "BLTree.h"
#include "Problems.h"
#include "GORoutines.h"
#include "IntegerVector.h"
#include "InputOutput.h"

/*----------------------------------------------------------------------------*/
void CheckBox (const PBLTREE & pbltwork, const PBLTREE & pbltfinal, 
		PBOX & pb, INTERVAL_VECTOR & GX, PBOX & pcF,
                HTscalar_FctPtr f, REAL & BetterUbf, VECTOR & BetterPoint,
		REAL Epsilon,
               const INTERVAL_VECTOR & XLim, INT WWidth, bool Draw,
	       bool MTest, bool EFatherGX, int & NEvalFx, int & NEvalFX, int & NEvalG, 
	       int & NDelMPT, int & NDelMT, int & NDelCTT, int & NDelFathGX);

void DivideBox(const PBLTREE & pbltwork, const PBLTREE & pbltfinal, 
		PBOX & pb, INTERVAL_VECTOR & GX, PBOX & pc,
               HTscalar_FctPtr f, const int & NDimDiv, 
               const int & NFrac, REAL & BetterUbf, VECTOR & BetterPoint,
               REAL Epsilon,
               const INTERVAL_VECTOR & XLim, INT WWidth, 
               bool Draw, bool MTest, bool EFatherGX,
	       int & NEvalFx, int & NEvalFX, int & NEvalG, 
	       int & NDelMPT, int & NDelMT, int & NDelCTT, int & NDelFathGX);

#endif /*__PROCESSBOX__*/
