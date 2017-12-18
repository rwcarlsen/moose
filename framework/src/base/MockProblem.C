/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "MockProblem.h"
#include "NonlinearSystem.h"

template <>
InputParameters
validParams<MockProblem>()
{
  InputParameters params = validParams<FEProblem>();
  params.addParam<int>("n_residual", 1, "number of residual evaluations for the mock solve");
  params.addParam<int>("n_jacobian", 1, "number of jacobian evaluations for the mock solve");
  return params;
}

MockProblem::MockProblem(const InputParameters & parameters)
  : FEProblem(parameters),
    _n_residual(getParam<int>("n_residual")),
    _n_jacobian(getParam<int>("n_jacobian"))
{
}

void
MockProblem::solve()
{
  for (int i = 0; i < _n_residual; i++)
    computeResidual(*_nl->currentSolution(), _nl->RHS());
  for (int i = 0; i < _n_jacobian; i++)
    computeJacobian(*_nl->currentSolution(), *getNonlinearSystem().sys().matrix, Moose::KT_ALL);
}
