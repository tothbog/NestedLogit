#include <stdlib.h>
#include <values.h>
#include <iostream> // .h was
#include "List.h"

/*----------------------------------------------------------------------------*/
LISTNODE::LISTNODE(PBOX & pb)
{
 pbox=pb;
 pnext=NULL;
 pb=NULL;
}

/*----------------------------------------------------------------------------*/
LISTNODE::~LISTNODE()             //Otherwise only the pointer pbox is deleted.
{
 if ( pbox != NULL )
    {
     delete pbox;
     pbox=NULL;
    } 
}

/*----------------------------------------------------------------------------*/
ostream & operator << (ostream & os, const LISTNODE & ln )
{
 if (ln.pbox !=NULL)
    os << "ListNode.*pbox="  << endl << *(ln.pbox) << endl;
 else
    os << "ListNode.pbox=NULL" << endl;

 os << "ListNode.pnext=" << ln.pnext; 
 
 return os;
}

/*----------------------------------------------------------------------------*/
/******************************************************************************/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
LIST::LIST()
{
 pfirst=NULL;
 nElements=0;
}

/*----------------------------------------------------------------------------*/
LIST::~LIST()
{
 PLISTNODE ptemp;

 if ( isEmpty() )
    {
     cerr << "LIST::~LIST, Empty List" << endl;
     exit(1);
    } 
 else
     while (pfirst!=NULL)
	   {
        ptemp=pfirst;
	    pfirst=pfirst->pnext;
	    delete ptemp;
	   }
}

/*----------------------------------------------------------------------------*/
bool LIST::isEmpty()
{
 if ( (pfirst==NULL && nElements!=0) || (pfirst!=NULL && nElements==0))
    {
     cerr << "LIST::isEmpty, nElements and pfirst do not match" << endl;
     exit(1);
    }
 return (pfirst==NULL);
}

/*----------------------------------------------------------------------------*/
int LIST::Dimension()
{
 return nElements;
}

/*----------------------------------------------------------------------------*/
void LIST::Insert(PBOX & pb)
{
 PLISTNODE pnew,plprev,plnext,pltemp;

 pnew = new LISTNODE(pb);  //I can not acces directly to pb because it is NULL.

 if (isEmpty())
    pfirst=pnew;
 else
    if ( Sup(pnew->pbox->FX)  <  Sup(pfirst->pbox->FX) )
       {
        pnew->pnext=pfirst;
        pfirst=pnew;
       }
    else
       {
        plprev=pfirst;
        plnext=plprev->pnext;
        while ( plnext!=NULL &&
       	       Sup(plnext->pbox->FX) <= Sup(pnew->pbox->FX) )
       	    {
       	       plprev=plnext;
       	       plnext=plnext->pnext;
       	      }
       	if (plnext!=NULL)
       	   {
       	    pltemp=plprev->pnext;
       	    plprev->pnext=pnew;
       	    pnew->pnext=pltemp;
       	   }
       	else
       	   plprev->pnext=pnew;
       }
 nElements++;
}

/*----------------------------------------------------------------------------*/
void LIST::GetFirst(PBOX & pb)
{
 if (isEmpty())
    {
     cerr << "LIST::GetFirst(), Empty List" << endl;
     exit(1);
    }
 if ( pfirst->pbox==NULL)
    {
     cerr << "LIST::GetFirst(), NULL Box" << endl;
     exit(1);     
    }    
 pb=pfirst->pbox;
 pfirst->pbox=NULL;
}


/*----------------------------------------------------------------------------*/
void LIST::RemoveFirst()
{
 PLISTNODE ptemp;
 
 if (isEmpty())
     {
      cerr << "LIST:RemoveFirst(), Empty List" << endl;
      exit(1);
     }
 ptemp=pfirst;
 pfirst=pfirst->pnext;
 delete ptemp; 														//Necesary.
 nElements--;
}

/*----------------------------------------------------------------------------*/
void LIST::Draw (const INTERVAL_VECTOR & XLim, INT WWidth, PCHAR color)
{
 PLISTNODE plcurr,plnext;
 
 if (isEmpty())
    {
     cerr << "LIST:Draw(), Empty List" << endl;
     exit(1);     
    } 
 else
    {
     plcurr=pfirst;
     plnext=plcurr->pnext;
     plcurr->pbox->Draw(XLim,WWidth,color);
     while (plnext!=NULL)
           {
            plcurr=plnext;
            plcurr->pbox->Draw(XLim,WWidth,color);
            plnext=plnext->pnext;
           }
    }   
}  									   //plcurr is deleted but not its contents.

/*----------------------------------------------------------------------------*/
void LIST::UnDraw ()
{
 PLISTNODE plcurr,plnext;
 
 if (isEmpty())
    {
     cerr << "LIST:Draw(), Empty List" << endl;
     exit(1);     
    } 
 else
    {
     plcurr=pfirst;
     plnext=plcurr->pnext;
     plcurr->pbox->UnDraw();
     while (plnext!=NULL)
           {
            plcurr=plnext;
            plcurr->pbox->UnDraw();
            plnext=plnext->pnext;
           }
    }   
}  									   //plcurr is deleted but not its contents.


/*----------------------------------------------------------------------------*/
REAL LIST::MinlbF()
{
 REAL MinlbF=DBL_MAX;
 
 if (isEmpty())
    {
     cerr << "LIST:MinlbF: Empty List" << endl;
     exit(1);     
    } 
 else
    {
     PLISTNODE plcurr,plnext;
     
     plcurr=pfirst;
     plnext=plcurr->pnext;
     if ( MinlbF > Inf(plcurr->pbox->FX) )
        MinlbF=Inf(plcurr->pbox->FX);
     while ( plnext != NULL )
           {
            plcurr=plnext;
	        if (MinlbF > Inf(plcurr->pbox->FX) )
        	   MinlbF=Inf(plcurr->pbox->FX);
             plnext=plnext->pnext;
            }
    }
 return MinlbF;
}


/*----------------------------------------------------------------------------*/
ostream & operator << (ostream & os, const LIST & l)
{
 PLISTNODE ptemp;

 os << "List.Dimension  = " << l.nElements << endl;
 os << "List ->" ;

 ptemp=l.pfirst;
 while (ptemp!=NULL)
       {
        os << endl << *ptemp ;
	    os << "->";
	    ptemp=ptemp->pnext;
       }
 os << "NULL"<< endl;

 return os;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
