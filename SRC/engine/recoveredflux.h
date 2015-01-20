// -*- C++ -*-
// $RCSfile: recoveredflux.h,v $
// $Revision: 1.5.10.1 $
// $Author: langer $
// $Date: 2013/11/08 20:44:45 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov. 
 */

#include <oofconfig.h>

#ifndef RECOVEREDFLUX_H
#define RECOVEREDFLUX_H

#include "common/doublevec.h"


//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//=\\=//

class RecoveredFlux {
private:
  DoubleVec *fvalue;
  int nfv;
public:
  RecoveredFlux(DoubleVec*);
  ~RecoveredFlux() { delete fvalue; }
  DoubleVec *get_flux_value() { return fvalue; };
  void average(DoubleVec*);
};

#endif // RECOVEREDFLUX
