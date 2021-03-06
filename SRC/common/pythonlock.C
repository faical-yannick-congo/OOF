// -*- C++ -*-
// $RCSfile: pythonlock.C,v $
// $Revision: 1.2.2.1 $
// $Author: langer $
// $Date: 2013/01/28 16:58:08 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov. 
 */

#include <oofconfig.h>

#include "common/pythonlock.h"

PyGILState_STATE acquirePyLock() {
  return PyGILState_Ensure();
}

void releasePyLock(PyGILState_STATE state) {
  PyGILState_Release(state);
}
