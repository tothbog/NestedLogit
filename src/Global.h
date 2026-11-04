#ifndef __GLOBAL__
//----------------------------------------------------------------------------
// VMGlobal.H
//
// (C) 1993 University of Karlsruhe
//
//----------------------------------------------------------------------------
#define __GLOBAL__

//----------------------------------------------------------------------------
// Klasse "index"
//
// Diese Dummy-Klasse hat einzig und allein den Zweck, den Konstruktor mit
// einem Parameter fuer rvector zu verschleiern; damit also keine direkte
// Typumwandlung von int in rvector moeglich ist. Dies ergaebe sonst in
// Anwenderprogrammen jede Menge Doppeldeutigkeiten, z. B.:
//  real a;
//  a=2*a; // ist dies nun:  rvector(int) * real
//                    oder:  real(double(int)) * real ???
//  So ist auf jeden Fall eine Umwandlungsstufe dazwischen geschoben!
//

#define index CXSC_index

class index {
  INT ind;
public:
  index() {}
  index(INT i) { ind=i; }
  inline INT _int() { return ind; }
};

#endif
