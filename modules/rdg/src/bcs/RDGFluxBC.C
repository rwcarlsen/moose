/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGFluxBC.h"

template<>
InputParameters validParams<RDGFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription("TODO.");
  params.addParam<unsigned int>("component", 0, "The componnet of the flux to use");
  params.addRequiredParam<UserObjectName>("flux", "Name of the boundary flux object to use");
  return params;
}

RDGFluxBC::RDGFluxBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _flux(getUserObject<BoundaryFluxBase>("flux"))
{
}

Real
RDGFluxBC::computeQpResidual()
{
  return _flux.getFlux()[_component] * _test[_i][_qp];
}

Real
RDGFluxBC::computeQpJacobian()
{
  return _flux.getJacobian()(_component, _component) * _phi[_j][_qp] * _test[_i][_qp];
}
