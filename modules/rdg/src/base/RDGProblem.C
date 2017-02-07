/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGProblem.h"
#include "RDGSystem.h"

template<>
InputParameters validParams<RDGProblem>()
{
  InputParameters params = validParams<FEProblemBase>();
  return params;
}

RDGProblem::RDGProblem(const InputParameters & parameters)
  : FEProblemBase(parameters)
{
  _nl = new RDGSystem(*this, "rdg0");
  _aux = new AuxiliarySystem(*this, "aux0");

  newAssemblyArray(*_nl);

  initNullSpaceVectors(parameters, *_nl);
}

RDGProblem::~RDGProblem()
{
  FEProblemBase::deleteAssemblyArray();

  delete _nl;
  delete _aux;
}

void
RDGProblem::solve()
{
  Moose::perf_log.push("solve()", "Execution");

  if (_solve)
    _nl->solve();

  if (_solve)
    _nl->update();

  Moose::perf_log.pop("solve()", "Execution");
}
