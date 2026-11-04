#include <iostream>
#include "Box.h"

extern int n;            // number of demand points

static INT BoxCont = 0;  // counter used to give each box an identifier

/*----------------------------------------------------------------------------*/
BOX::BOX (INT nn)
{
 if (nn<0)
    {
     cerr << "BOX::BOX(int) negative size" << endl;
     exit(1);
    }
 Resize(X,nn);
 M=false;
 BoxCont++;
 NBox=BoxCont;
 FatherNBox = -1;
 feas=false;
 Resize(fe, 2*n); 
 Clear(fe);
}

/*----------------------------------------------------------------------------*/
BOX::BOX (const INTERVAL_VECTOR & IV)
{
 int NDim = Dimension(IV);
 Resize(X,NDim);
 X=IV;
 M=false;
 BoxCont++;
 NBox=BoxCont; 
 FatherNBox = -1;
 feas=false;
 Resize(fe, 2*n); Clear(fe);
}

/*----------------------------------------------------------------------------*/
BOX::BOX(const BOX & B)
{
 int NDim = Dimension(B.X);
 Resize(X,NDim);
 X=B.X;
 FX=B.FX;
 M=B.M;
 NBox=B.NBox;
 FatherNBox=B.FatherNBox;
 feas=B.feas;
 fe=B.fe;
}

/*----------------------------------------------------------------------------*/
BOX::BOX(bool son, const BOX & B)
{
 X=B.X;
 FX=B.FX;
 M=B.M;
 if (son)
 {
   BoxCont++;
   NBox=BoxCont;
   FatherNBox = B.NBox;
 }
 else
 {
   NBox=B.NBox;
   FatherNBox = B.FatherNBox;
 }  
 feas=B.feas;
 fe=B.fe;
}

/*----------------------------------------------------------------------------*/
BOX::~BOX ()
{
 Resize(X,0);
 Resize(fe,0);
}

/*----------------------------------------------------------------------------*/
int Dimension (const BOX & B)
{
 return Dimension(B.X);
}

/*----------------------------------------------------------------------------*/
REAL EvalpF(const BOX & B, int Index, REAL BetterUbf)
{ 
 return Sup ( ( Hull(BetterUbf) - Hull(Inf (B.FX)) ) / 
              Hull(Diam(B.FX)) );
}

/*----------------------------------------------------------------------------*/
/* Write the (x1,x2) projection of the box, scaled to a WWidth x WWidth      */
/* window (graphical output).                                                  */
/*----------------------------------------------------------------------------*/
void BOX::Draw (const INTERVAL_VECTOR & XLim, INT WWidth, PCHAR color, ostream & out )
{
 REAL FYDiam,FXDiam;
 
 FXDiam = Diam(XLim(1));
 FYDiam = Diam(XLim(2));
 
 out << "PintaBox" << endl;
 out << (Inf(X(1))-Inf(XLim(1)))  / FXDiam * (REAL)(WWidth-1) << endl; 
 out << ((Inf(X(2))-Inf(XLim(2))) / FYDiam * (REAL)(-1.0) + (REAL)1.0 ) *
         (REAL)(WWidth-1) << endl;
 out << (Sup(X(1))-Inf(XLim(1)))  / FXDiam * (REAL)(WWidth-1) << endl; 
 out << ((Sup(X(2))-Inf(XLim(2))) /FYDiam  * (REAL)(-1.0) + (REAL)1.0 ) *
         (REAL)(WWidth-1) << endl; 
 out << color << endl;
 out << NBox << endl;
}

/*----------------------------------------------------------------------------*/
void BOX::Draw (const INTERVAL_VECTOR & XLim, INT WWidth, PCHAR color )
{
 Draw(XLim, WWidth, color, cout);
}

/*---------------------------------------------------------------------------*/
void BOX::UnDraw ()
{
 cout << "BorraBox" << endl;
 cout << NBox << endl;
}

/*----------------------------------------------------------------------------*/
ostream & operator << (ostream & os, const BOX & B)
{
  os << "\tX           = " << B.X  <<endl;
  os << "\tw(X)        = " << MaxDiam(B.X)<< endl;
  os << "\tFX          = " << B.FX << endl;
  os << "\tMonot       = " << B.M << endl;
  os << "\tNBox        = " << B.NBox << endl;
  os << "\tFatherNBox  = " << B.FatherNBox << endl;
  os << "\tfeasible    = " << B.feas << endl; 
  os << B.fe;
  os << endl;
  return os;
}

/*---------------------------------------------------------------------------*/
ostream & operator << (ostream & os, PBOX & pb)
{
 if (pb==NULL)
   os << "NULL" << endl;
 else
   os << *pb;
 
 return os;
}
