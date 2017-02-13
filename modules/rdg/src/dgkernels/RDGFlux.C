/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGFlux.h"
#include "InternalSideFluxBase.h"

template<>
InputParameters validParams<RDGFlux>()
{
  InputParameters params = validParams<DGKernel>();
  params.addClassDescription("A dgkernel for the advection equation using a cell-centered finite volume method.");
  params.addParam<unsigned int>("component", 0, "Which component of the flux vector to use");
  params.addRequiredParam<UserObjectName>("flux", "Name of the internal side flux object to use");
  return params;
}

RDGFlux::RDGFlux(const InputParameters & parameters) :
    DGKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _flux(getUserObject<InternalSideFluxBase>("flux"))
{
}

RDGFlux::~RDGFlux()
{
}

Real
RDGFlux::computeQpResidual(Moose::DGResidualType type)
{
  const auto & flux = _flux.getFlux();

  // std::cerr << "flux = " << flux[_component] << ", " << _test[_i][_qp] << std::endl;

  // distribute the contribution to the current and neighbor elements
  switch (type)
  {
    case Moose::Element:
      return flux[_component] * _test[_i][_qp];

    case Moose::Neighbor:
      return -flux[_component] * _test_neighbor[_i][_qp];
  }

  return 0.0;
}

Real
RDGFlux::computeQpJacobian(Moose::DGJacobianType type)
{
  const auto & fjac1 = _flux.getJacobian(Moose::Element);
  const auto & fjac2 = _flux.getJacobian(Moose::Neighbor);

  // distribute the contribution to the current and neighbor elements
  switch (type)
  {
    case Moose::ElementElement:
      return fjac1(_component, _component) * _phi[_j][_qp] * _test[_i][_qp];

    case Moose::ElementNeighbor:
      return fjac2(_component, _component) * _phi_neighbor[_j][_qp] * _test[_i][_qp];

    case Moose::NeighborElement:
      return -fjac1(_component, _component) * _phi[_j][_qp] * _test_neighbor[_i][_qp];

    case Moose::NeighborNeighbor:
      return -fjac2(_component, _component) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
  }

  return 0.0;
}
