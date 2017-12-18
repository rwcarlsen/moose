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

#ifndef MOCKPROBLEM_H
#define MOCKPROBLEM_H

#include "FEProblem.h"

class FEProblem;

template <>
InputParameters validParams<FEProblem>();

class MockProblem : public FEProblem
{
public:
  MockProblem(const InputParameters & parameters);

  virtual bool converged() override { return true; }

  virtual void solve() override;

  virtual ~MockProblem(){};

private:
  int _n_residual = 0;
  int _n_jacobian = 0;
};

#endif /* MOCKPROBLEM_H */
