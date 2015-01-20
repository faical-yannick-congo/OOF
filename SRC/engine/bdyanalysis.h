// -*- C++ -*-
// $RCSfile: bdyanalysis.h,v $
// $Revision: 1.8 $
// $Author: langer $
// $Date: 2014/09/27 21:40:40 $

/* This software was produced by NIST, an agency of the U.S. government,
 * and by statute is not subject to copyright in the United States.
 * Recipients of this software assume all responsibilities associated
 * with its operation, modification and maintenance. However, to
 * facilitate maintenance we ask that before distributing modified
 * versions of this software, you first contact the authors at
 * oof_manager@nist.gov.
 */

#include <oofconfig.h>

#ifndef BDYANALYSIS_H
#define BDYANALYSIS_H

#include <vector>

class EdgeSet;
class FEMesh;
class Flux;
class OutputValue;

OutputValue integrateFlux(const FEMesh*, const Flux*, const EdgeSet*);
OutputValue averageField(const FEMesh*, const Field*, const EdgeSet*);

#endif // BDYANALYSIS_H
