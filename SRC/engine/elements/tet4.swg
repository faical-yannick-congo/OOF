// -*- C++ -*-
// $RCSfile: tet4.swg,v $
// $Revision: 1.2.10.1 $
// $Author: langer $
// $Date: 2014/09/27 22:34:28 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov. 
 */

#ifndef TET4_SWIG
#define TET4_SWIG

%module tet4
%include "engine/typemaps.swg"
%{
  extern void tet4init();
%}
void tet4init();

%pragma(python) include="tet4.spy"

#endif // TET4_SWIG
