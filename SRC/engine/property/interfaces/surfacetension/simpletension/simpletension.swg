// -*- C++ -*-
// $RCSfile: simpletension.swg,v $
// $Revision: 1.2.16.1 $
// $Author: langer $
// $Date: 2013/11/08 20:45:57 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov. 
 */

%module simpletension

%{
#include "engine/property/interfaces/surfacetension/simpletension/simpletension.h"
%}
%include "common/typemaps.swg"
%include "engine/propertyhdr.swg"

class SimpleTension : public Property
{
public:
  SimpleTension(PyObject *self, char *, double, double);
};

%pragma(python) include="simpletension.spy"
