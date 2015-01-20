// -*- C++ -*-
// $RCSfile: general_nonlinear_elasticity.swg,v $
// $Revision: 1.4 $
// $Author: gdogan $
// $Date: 2010/12/06 17:31:56 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov.
 */

%module nonlinear

%{
#include "engine/property/elasticity/nonlinear/general_nonlinear_elasticity.h"
%}

%include "common/typemaps.swg"
%include "engine/propertyhdr.swg"
%pragma(python) include="general_nonlinear_elasticity.spy"

class GeneralNonlinearElasticityNoDeriv : public Property {
};

class GeneralNonlinearElasticity : public GeneralNonlinearElasticityNoDeriv {
};

class TestGeneralNonlinearElasticityNoDeriv : public GeneralNonlinearElasticityNoDeriv {
public:
  TestGeneralNonlinearElasticityNoDeriv(PyObject *reg, char *name, int testno);
};

class TestGeneralNonlinearElasticity : public GeneralNonlinearElasticity {
public:
  TestGeneralNonlinearElasticity(PyObject *reg, char *name, int testno);
};