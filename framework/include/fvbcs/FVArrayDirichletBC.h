//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVBoundaryCondition.h"
#include "MooseVariableInterface.h"

/**
 * Base class for
 */
class FVArrayDirichletBC : public FVBoundaryCondition,
                           public MooseVariableInterface<RealEigenVector>
{
public:
  /**
   * Class constructor.
   * @param parameters The InputParameters for the object
   * @param nodal Whether this BC is applied to nodes or not
   */
  FVArrayDirichletBC(const InputParameters & parameters);

  static InputParameters validParams();

  virtual RealEigenVector boundaryValue(const FaceInfo & fi) const;

private:
  const RealEigenVector & _val;
};
