// -*- C++ -*-
// $RCSfile: orientmapproperty.swg,v $
// $Revision: 1.1.26.1 $
// $Author: langer $
// $Date: 2013/11/08 20:46:07 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov. 
 */

#ifndef ORIENTMAPPROPERTY_SWG
#define ORIENTMAPPROPERTY_SWG

%module orientationmapproperty
%include "common/typemaps.swg"
%include "engine/propertyhdr.swg"
%{
#include "orientationmap/orientmapproperty.h"
%}

%pragma(python) include="orientationmap/orientmapproperty.spy"

class OrientationMapProp : public Property {
public:
  OrientationMapProp(PyObject *registry, char *name);
};

#endif // ORIENTMAPPROPERTY_SWG
