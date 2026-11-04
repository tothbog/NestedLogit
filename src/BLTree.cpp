#include <stdlib.h>
#include <values.h>
#include <iostream>
#include "BLTree.h"
#include "Problems.h"

/*----------------------------------------------------------------------------*/
BLTREENODE::BLTREENODE(PBOX & pb, REAL pbKey)
{
 pList = new(LIST);
 pList->Insert(pb);
 Key=pbKey;
 Balance=EQUAL;
 pleft=NULL;
 pright=NULL;
}	

/*----------------------------------------------------------------------------*/
BLTREENODE::~BLTREENODE()
{
 delete pList;
 pList=NULL;
}

/*----------------------------------------------------------------------------*/
/*Swap the List and Key of the nodes, but the structure of the BLTree does not*/
/*change.																	  */
/*----------------------------------------------------------------------------*/                                                                    
void SwapData(const PBLTREENODE & node1, const PBLTREENODE & node2)
{
 PLIST  pListTemp;
 REAL  KeyTemp;
 
 pListTemp       = node1->pList;
 KeyTemp    	= node1->Key;
 node1->pList 	= node2->pList;
 node1->Key  	= node2->Key;
 node2->pList 	= pListTemp;
 node2->Key  	= KeyTemp;
}

/*----------------------------------------------------------------------------*/
/*Swap the data and copy the Balance and pointers of the BLTreeNodes.         */
/*Node2 will be removed with data of node1 after exchange.                    */
/*----------------------------------------------------------------------------*/                                                                    
void Exchange(const PBLTREENODE & node1, const PBLTREENODE & node2)
{
 SwapData(node1,node2);
 
 node1->Balance	= node2->Balance;
 node1->pleft	= node2->pleft;
 node1->pright	= node2->pright;

 node2->pleft	= NULL;
 node2->pright	= NULL;
}

/*----------------------------------------------------------------------------*/
ostream & operator << (ostream & os, const BLTREENODE & bltn)
{
 os << "BLTreeNode.Key\t\t= " 	<< bltn.Key               << endl;
 os << "BLTreeNode.Balance\t= " << (int) bltn.Balance     << endl;
 os << "BLTreeNode.*pList\t= " 	<<  endl << *(bltn.pList) << endl;
 os << "BLTreeNode.pleft\t= " 	<< bltn.pleft             << endl; 
 os << "BLTreeNode.pright\t= " 	<< bltn.pright            << endl;  
 
 return os;
}

/*----------------------------------------------------------------------------*/
/******************************************************************************/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
BLTREE::BLTREE()
{
 root=NULL;
 nElements=0;
 MaxnElements=0;
 KeyCut=-DBL_MAX; 
}

/*----------------------------------------------------------------------------*/
void DeleteBLSubTree(PBLTREENODE & root)
{
 if (root!=NULL)
 	{
 	 DeleteBLSubTree(root->pleft);
 	 DeleteBLSubTree(root->pright);
 	 delete root;
 	 root=NULL;
 	}
}

/*----------------------------------------------------------------------------*/
BLTREE::~BLTREE()
{
 DeleteBLSubTree(root);
}

/*----------------------------------------------------------------------------*/
bool BLTREE::isEmpty()
{
 if ( (root==NULL && nElements!=0) || (root!=NULL && nElements==0))
    {
     cerr << "BLTree::isEmpty: nElements and root do not match" << endl;
     exit(1);
    }
 return (root==NULL);
}

/*----------------------------------------------------------------------------*/
int BLTREE::Dimension()
{
 return nElements;
}

/*----------------------------------------------------------------------------*/
int BLTREE::MaxDimension()
{
 return MaxnElements;
}

/*----------------------------------------------------------------------------*/
REAL BLTREE::GetKeyCut()
{
 return KeyCut;
}

/*----------------------------------------------------------------------------*/
/*Rotate left a binary Tree                                                   */
/*----------------------------------------------------------------------------*/
void RotateLeft(PBLTREENODE & root)
{
 PBLTREENODE Temp;

 if (root==NULL)
    {
     cerr << "BLTREE::RotateLeft :";
     cerr << " Can not rotate to left an empty tree " << endl;
     exit(1);
    }
 else
    if (root->pright==NULL)
       {
        cerr << "BLTREE::RotateLeft : Can not rotate to left" << endl;
        exit(1);
       }
    else
       {
        Temp         = root->pright;
        root->pright = Temp->pleft;
        Temp->pleft  = root;
        root=Temp;
       }
}

/*----------------------------------------------------------------------------*/
/*Rotate right a binary Tree                                                  */
/*----------------------------------------------------------------------------*/
void RotateRight(PBLTREENODE & root)
{
 PBLTREENODE Temp;
 
 if (root==NULL)
    {
     cerr << "BLTREE::RotateRight :";
     cerr << " Can not rotate to right an empty tree " << endl;
     exit(1);
    }
 else
    if (root->pleft==NULL)
       {
        cerr << "BLTREE::RotateRight : Can not rotate to right" << endl;
        exit(1);
       }
    else
       {
        Temp         = root->pleft;
        root->pleft  = Temp->pright;
        Temp->pright = root;
        root=Temp;
       }
}

/*----------------------------------------------------------------------------*/
/*Insert in the left branch maintaining the global balance. Taller noticed if */
/*the insertion has increased the heigh or not, so the balance factors can be */
/*changed appropiately and rotations may be done. Taller has efects only as   */
/* return value. We initiate it to false.                                     */
/*----------------------------------------------------------------------------*/
void LeftInsertBalance(PBLTREENODE & root, bool & taller)
{
 PBLTREENODE Left, LeftRight;

 Left=root->pleft;
 switch (Left->Balance)
        {
         case LEFT:
              root->Balance=EQUAL;
              Left->Balance=EQUAL;
              RotateRight(root);
              taller=false;
              break;
         case EQUAL:
              cerr << "BLTREE:LeftInsertBalance : The tree is already balanced";
              cerr << endl;
              exit(1);
              break;
         case RIGHT:
              LeftRight=Left->pright;
              switch (LeftRight->Balance)
                     {
                      case LEFT:
                           root->Balance=RIGHT;
                           Left->Balance=EQUAL;
                           break;
                      case EQUAL:
                           root->Balance=EQUAL;
                           Left->Balance=EQUAL;
                           break;
                      case RIGHT:
                           root->Balance=EQUAL;
                           Left->Balance=LEFT;
                           break;
                     }
              LeftRight->Balance=EQUAL;
              RotateLeft(root->pleft);
              RotateRight(root);
              taller=false;
              break;
        }
}

/*---------------------------------------------------------------------------*/
/*Analogous to previous one but on the right.                                */
/*---------------------------------------------------------------------------*/
void RightInsertBalance(PBLTREENODE & root, bool & taller)
{
 PBLTREENODE Right, RightLeft;

 Right=root->pright;
 switch (Right->Balance)
        {
         case LEFT:
              RightLeft=Right->pleft;
              switch (RightLeft->Balance)
                     {
                      case LEFT:
                           root->Balance=EQUAL;
                           Right->Balance=RIGHT;
                           break;
                      case EQUAL:
                           root->Balance=EQUAL;
                           Right->Balance=EQUAL;
                           break;
                      case RIGHT:
                           root->Balance=LEFT;
                           Right->Balance=EQUAL;
                           break;
                     }
              RightLeft->Balance=EQUAL;
              RotateRight(root->pright);
              RotateLeft(root);
              taller=false;
              break;
         case EQUAL:
              cerr << "BLTREE:RightInsertBalance: The tree is already balanced";
              cerr << endl;
              exit(1);
              break;
         case RIGHT:
              root->Balance=EQUAL;
              Right->Balance=EQUAL;
              RotateLeft(root);
              taller=false;
              break;
        }
}

/*---------------------------------------------------------------------------*/
/*Do the actual insertion in the BLTREE.                                     */
/*---------------------------------------------------------------------------*/ 
void Insert1 (PBLTREENODE & root, PBOX & pb, REAL pbKey, bool & taller)
{
 if (root==NULL)
    {
     root = new BLTREENODE(pb,pbKey);
     taller=true;
    }
 else
     if (pbKey==root->Key)
        root->pList->Insert(pb);
     else
         if (pbKey < root->Key)
            {
             Insert1(root->pleft, pb, pbKey, taller);
             if (taller)
                switch (root->Balance)
                       {
                        case LEFT:
                                 LeftInsertBalance(root,taller);
                                 break;
                        case EQUAL:
                                 root->Balance=LEFT;
                                 break;
                        case RIGHT:
                                 root->Balance=EQUAL;
                                 taller=false;
                                 break;
                       }
            }
         else
            {
             Insert1(root->pright, pb, pbKey, taller);
             if (taller)
                switch (root->Balance)
                       {
                        case LEFT:
                                 root->Balance=EQUAL;
                                 taller=false;
                                 break;
                        case EQUAL:
                                 root->Balance=RIGHT;
                                 break;
                        case RIGHT:
                                 RightInsertBalance(root,taller);
                                 break;
                       }
            }
}

/*---------------------------------------------------------------------------*/
void BLTREE::Insert(PBOX & pb, REAL pbKey)
{
 bool taller=false;

 if (Inf(pb->FX) > KeyCut)
    KeyCut=Inf(pb->FX); 
    	
 nElements++;

 if (nElements > MaxnElements)
	 MaxnElements=nElements;

 Insert1(root,pb,pbKey,taller);
}


/*---------------------------------------------------------------------------*/
/*Update the balance factors and do the necesary rotations when a left node  */
/*is removed.                                                                */
/*---------------------------------------------------------------------------*/
void LeftDelBalance (PBLTREENODE & root, bool & shorter)
{
 PBLTREENODE Left, LeftRight;

 Left=root->pleft;
 switch (Left->Balance)
        {
         case LEFT:
              root->Balance=Left->Balance=EQUAL;
              RotateRight(root);
              break;
         case EQUAL:
              root->Balance=LEFT;
              Left->Balance=RIGHT;
              RotateRight(root);
              shorter=false;
              break;
         case RIGHT:
              LeftRight=Left->pright;
              switch (LeftRight->Balance)
                     {
                      case LEFT:
                           root->Balance=RIGHT;
                           Left->Balance=EQUAL;
                           break;
                      case EQUAL:
                           root->Balance=Left->Balance=EQUAL;
                           break;
                      case RIGHT:
                           root->Balance=EQUAL;
                           Left->Balance=LEFT;
                           break;
                     }
              LeftRight->Balance=EQUAL;
              RotateLeft(root->pleft);
              RotateRight(root);
              break;
        }
}

/*---------------------------------------------------------------------------*/
/*Analogous to previous one but on the right.                                */
/*---------------------------------------------------------------------------*/
void RightDelBalance (PBLTREENODE & root, bool & shorter)
{
 PBLTREENODE Right, RightLeft;

 Right=root->pright;
 switch (Right->Balance)
        {
         case LEFT:
              RightLeft=Right->pleft;
              switch (RightLeft->Balance)
                     {
                      case LEFT:
                           root->Balance=EQUAL;
                           Right->Balance=RIGHT;
                           break;
                      case EQUAL:
                           root->Balance=Right->Balance=EQUAL;
                           break;
                      case RIGHT:
                           root->Balance=LEFT;
                           Right->Balance=EQUAL;
                           break;
                     }
              RightLeft->Balance=EQUAL;
              RotateRight(root->pright);
              RotateLeft(root);
              break;
         case EQUAL:
              root->Balance=RIGHT;
              Right->Balance=LEFT;
              RotateLeft(root);
              shorter=false;
              break;
         case RIGHT:
              root->Balance=Right->Balance=EQUAL;
              RotateLeft(root);
              break;
        }
}

/*---------------------------------------------------------------------------*/
/*Delete the ListNode with Key=pbKey. When Pop=False, the whole list is      */
/*removed, otherwise only the first list element is removed. The variable    */
/*shorter is initially True. When it becomes False, no further changes are   */
/*needed. nDeleted return the number of deleted boxes.                                                                   */
/*---------------------------------------------------------------------------*/ 
void Del1(PBLTREENODE & root, REAL pbKey, bool & shorter, const bool & Pop, 
int & nDeleted)
{
 PBLTREENODE Temp;
 INTERVAL kk=Hull(0.0);

 if (root==NULL)
    {
     cerr << "DelNode: Empty BLTREE" << endl;
     cerr << 10.0/kk << endl;
    }     
 else       
    if (pbKey==root->Key)
       if (Pop && root->pList->Dimension()>1) /*Otherwise the list is removed*/
          {
           nDeleted=1;
           root->pList->RemoveFirst();
	       shorter=false;
          } 
       else
          if (root->pleft==NULL && root->pright==NULL)
             {
              nDeleted=root->pList->Dimension();
              delete root;
              root=NULL;
              shorter=true;
             } 
          else
             if (root->pleft==NULL)
                {
                 Temp=root->pright;
                 Exchange(root,Temp); 
                 shorter=true;
                 nDeleted=Temp->pList->Dimension();
                 delete Temp;
                 Temp=NULL;
                } 
             else
                if (root->pright==NULL)
                   {
                    Temp=root->pleft;
                    Exchange(root,Temp);
                    shorter=true;
                    nDeleted=Temp->pList->Dimension();
                    delete Temp;
                    Temp=NULL;
                   } 
                else
                   {
                    for (Temp=root->pright;Temp->pleft!=NULL; Temp=Temp->pleft);
                    SwapData(root,Temp);  
                    Del1(root->pright, pbKey, shorter, Pop, nDeleted);
                    if (shorter)
                       switch (root->Balance)
                              {
                               case LEFT:
                                    LeftDelBalance(root, shorter);
                                    break;
                               case EQUAL:
                                    root->Balance=LEFT;
                                    shorter=false;
                                    break;
                               case RIGHT:
                                    root->Balance=EQUAL;
                                    break;
                              }
                   } 
     else
        if (pbKey<root->Key)
           {
            Del1(root->pleft, pbKey, shorter, Pop, nDeleted);
            if (shorter)
               switch (root->Balance)
                      {
                       case LEFT:
                            root->Balance=EQUAL;
                            break;
                       case EQUAL:
                            root->Balance=RIGHT;
                            shorter=false;
                            break;
                       case RIGHT:
                            RightDelBalance(root, shorter);
                            break;
                      }              
           }
        else
           {
            Del1(root->pright, pbKey, shorter, Pop, nDeleted);
            if (shorter)
                switch (root->Balance)
                      {
                       case LEFT:
                            LeftDelBalance(root, shorter);
                            break;
                       case EQUAL:
                            root->Balance=LEFT;
                            shorter=false;
                            break;
                       case RIGHT:
                            root->Balance=EQUAL;
                            break;
                      }
           }
}

/*---------------------------------------------------------------------------*/ 
void BLTREE::Del(REAL pbKey, const bool & Pop)
{
 bool shorter=true;
 int nDeleted=0;

 Del1(root, pbKey, shorter, Pop, nDeleted);
 nElements-=nDeleted;
 if (isEmpty())
 	KeyCut=-DBL_MAX;
   	//KeyCut is updated in Insert, by the CutOffTest and here.
} 

/*---------------------------------------------------------------------------*/
bool SearchGTKeyBox (const PBLTREENODE & root, REAL & KeyCut, 
                     const REAL & KeyLimit, REAL & ReturnedKey,
                     const INTERVAL_VECTOR & XLim, INT WWidth, 
				     bool Draw, PCHAR color)
{
 PBLTREENODE Temp;

 Temp=root;
 while (Temp!=NULL && !(Temp->Key > KeyLimit))
       {
        KeyCut=Temp->Key;
        Temp=Temp->pright;
       }
 if (Temp==NULL)
    return false;

 ReturnedKey=Temp->Key;
 if (Draw)
    {
     Temp->pList->UnDraw();
     Temp->pList->Draw(XLim,WWidth,color);
    }
 return true;
}

/*---------------------------------------------------------------------------*/
void BLTREE::DelGTKey (const REAL & KeyLimit, const INTERVAL_VECTOR & XLim, 
                       INT WWidth, bool Draw, PCHAR color)
{
 REAL ReturnedKey;
 bool Pop=false;

 while (SearchGTKeyBox (root, KeyCut, KeyLimit, ReturnedKey,
                        XLim, WWidth, Draw, color) )
       Del(ReturnedKey, Pop);
}

/*---------------------------------------------------------------------------*/ 
void BLTREE::GetLwKeyBox(PBOX & pb)
{
 PBLTREENODE Temp;

 if (isEmpty())
 	{
 	 cerr << "BLTREE::GetLwKeyBox empty BLTree" << endl;
     exit(1);
 	}
 Temp=root;
 while (Temp->pleft!=NULL)
 	   Temp=Temp->pleft;
 Temp->pList->GetFirst(pb);
 Del(Temp->Key,true);
}


/*----------------------------------------------------------------------------*/
void MinlbFNodes ( PBLTREENODE root, REAL & MinlbF)
{
 if (root!=NULL)
    {
     REAL MinTemp;
     
     MinTemp=root->pList->MinlbF();
     if ( MinlbF > MinTemp )
        MinlbF=MinTemp;
     MinlbFNodes(root->pleft,MinlbF);
     MinlbFNodes(root->pright,MinlbF);     
    } 
}

/*----------------------------------------------------------------------------*/
REAL BLTREE::MinlbF()
{
 REAL MinlbF=DBL_MAX;

 MinlbFNodes(root,MinlbF);
 return MinlbF;
}



/*----------------------------------------------------------------------------*/
void PrintBLTreeNodes( PBLTREENODE root, ostream & os)
{
 if (root!=NULL)
    {
     os << *root << endl;
     PrintBLTreeNodes(root->pleft, os);
     PrintBLTreeNodes(root->pright, os);
    }
 else
 	os << "N" << endl;
}


/*----------------------------------------------------------------------------*/
ostream & operator << (ostream & os, const BLTREE & blt)
{
 os << "BLTree.Dimension\t= " << blt.nElements    << endl;
 os << "BLTree.MaxDimension\t= " << blt.MaxnElements << endl;
 os << "BLTree.KeyCut\t\t= " << blt.KeyCut       << endl;
 os << "--------------------BLTREE Nodes-------------------------" << endl;
 
 PrintBLTreeNodes(blt.root, os);
 return os;
}


/*---------------------------------------------------------------------------*/
/*****************************************************************************/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
BLTREEP_VECTOR::BLTREEP_VECTOR()
{
 nElements=0;
 theElements=NULL;
 MaxElements=0;
}

/*---------------------------------------------------------------------------*/
BLTREEP_VECTOR::BLTREEP_VECTOR (int n)
{
 int i;
 
 nElements   = n;
 MaxElements=0;
 theElements = new PBLTREE[n];
 CheckMemory (theElements);

 for (i=0;i<n;i++)
     {
 	  theElements[i] = new BLTREE;
 	  CheckMemory (theElements[i]);
 	 } 
}

/*---------------------------------------------------------------------------*/
BLTREEP_VECTOR::~BLTREEP_VECTOR()
{
 int i;

 for (i=0;i<nElements;i++)
     delete theElements[i];
 delete [] theElements;
}

/*---------------------------------------------------------------------------*/
int BLTREEP_VECTOR::Dimension()  
{
 return nElements;
}

/*---------------------------------------------------------------------------*/
PBLTREE & BLTREEP_VECTOR::operator () (int i)  
{
 return theElements[i-1];
}
/*---------------------------------------------------------------------------*/
int BLTREEP_VECTOR::MaxLength()  
{
 return MaxElements;
}

/*---------------------------------------------------------------------------*/
void BLTREEP_VECTOR::UpdateMaxLength()  
{
 int nodes=0, i, n = Dimension();
 for( i = 1; i <= n; i++ )
       nodes += theElements[i-1]->Dimension();
 if (nodes > MaxElements) MaxElements=nodes;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/ 
