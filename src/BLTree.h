/*
Purpose: Definition of a balanced binary tree where each node contains a Key
    and a linked list with the Boxes with the same Key. A vector of BLTREEs
    is also defined. It is used usually to keep together a working and final 
    data structures.

                 --------------------------------------------
Class BLTREENODE:
Stores:
      pList   : pointer to a List containing Boxes with the same Key.
      Key     : Defines the order in the BLTREE.
      Balance : Determines the balance of the tree which has this node
                as root.
      pleft   : Pointer to the left subtree.
      pright  : Pointer to the right subtree.

functions:

 Public:
  BLTREENODE(PBOX, REAL) : constructor. Insert a PBOX in the List and
  		initialises Key to the REAL. BLTREENODE() is not allowed.
  ~BLTREENDE : Delete de memory asigned to pList.	
  Assignament opperator overloaded by the compiler, i.e, all data is
  	copied (not the data pointed by the pointers).

                 --------------------------------------------
Class BLTREE: friend of BLTREENODE.
Stores:
      root : root of the tree.
      KeyCut : Max Inf(X.FX) for all X in BLTree.
      nElements : Number of nodes in the tree.
      MaxnElements : Max number of nElements.

functions:

 Public:
   BLTREE()		: Constructor.
   ~BLTREE()	: Destructor.
   isEmpty()	: To know if the BLTree is empty.
   Dimension()  : Return the number of Boxes in the BLTree.
   MaxDimension(): Return the max number of Boxes stored in the BLTree.
   GetKeyCut  (): Return the KeyCut value.    
   Insert (PBOX pb, REAL pbKey): Insert a Box in the BLTree based on pbKey.
   Del(REAL pbKey, bool Pop): Delete a LIST or the first ListNode of the 
   		BLTREENODE with Key=pbKey, depending on Pop. If the LIST becomes empty
   		the BLTREE node is removed and the BLTREE is balanced, when necessary.
   DelGTKey (const REAL & KeyLimit, const INTERVAL_VECTOR & XLim, 
                   INT WWidth, bool Draw, PCHAR color): Delete the BLTreeNodes
        with Key > KeyLimit. 
   GetLwKeyBox	(PBOX & pb): Obtain the box with smaller Key in the BLTREE and 
   		the LISTNODE containg it is is removed.
   MinlbF(): Return the minimum value of Inf(FX) in the BLTree.


                 --------------------------------------------
Class BLTREEP_VECTOR:
Stores:             
	nElements: Number of pointers to BLTREE stored in the vector.
	theElements: Vector of pointers to BLTREEs.

  public:
   BLTREEP_VECTOR()		: Constructor. No memory allocated.
   BLTREEP_VECTOR(int n): Construct a vector of n pointers to BLTREE and 
   		construct each BLTREE.
   ~BLTREEP_VECTOR()	: Destructor.
   PBLTREE & operator()(int i) : Changes the initial index from 0 to 1.
	
*******************************************************************************/
#ifndef __BLTREE__
#define __BLTREE__

#include "List.h"

typedef enum BALANCEFACTOR {EQUAL, LEFT, RIGHT} BALANCETYPE;

class BLTREENODE
{
 friend class BLTREE;

 public:
	PLIST        pList;
	REAL         Key;
	BALANCETYPE  Balance;
	BLTREENODE * pleft;
	BLTREENODE * pright;
	
 public:
 	BLTREENODE(PBOX & pb, REAL pbKey);
 	~BLTREENODE();
    friend ostream & operator << (ostream & os, const BLTREENODE & bltn);
 private:
 	BLTREENODE(); //Default constructor is not allowed.		
};

/*----------------------------------------------------------------------------*/
typedef BLTREENODE * PBLTREENODE;
typedef PBLTREENODE * PPBLTREENODE;

/*----------------------------------------------------------------------------*/	
class BLTREE
{
 private:
	PBLTREENODE root;
	REAL		KeyCut;
	INT         nElements;
	INT			MaxnElements;
	
 public:
	BLTREE();
	~BLTREE();
	bool isEmpty	();
	int Dimension	();
	int MaxDimension();
	REAL GetKeyCut  ();
	void Insert	(PBOX & pb, REAL pbKey);
	void Del 	(REAL pbKey, const bool & Pop);
	void DelGTKey   (const REAL & KeyLimit, const INTERVAL_VECTOR & XLim, 
                         INT WWidth, bool Draw, PCHAR color);
	void GetLwKeyBox(PBOX & pb);
	REAL MinlbF     ();
	friend ostream & operator << (ostream &, const BLTREE &);
};

typedef BLTREE * PBLTREE;
typedef PBLTREE * PPBLTREE;

/*----------------------------------------------------------------------------*/	
class BLTREEP_VECTOR
{
 private:
 	int     nElements;
 	PPBLTREE theElements;
	int	MaxElements;
 	
 public:
 	BLTREEP_VECTOR();
 	BLTREEP_VECTOR(int n);
 	~BLTREEP_VECTOR();
 	int Dimension();
 	PBLTREE & operator()(int i);
	int MaxLength();
	void UpdateMaxLength();
};


#endif /*__BLTREE__*/
