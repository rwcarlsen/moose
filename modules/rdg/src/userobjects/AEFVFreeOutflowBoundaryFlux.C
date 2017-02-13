/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "AEFVFreeOutflowBoundaryFlux.h"

template<>
InputParameters validParams<AEFVFreeOutflowBoundaryFlux>()
{
  InputParameters params = validParams<BoundaryFluxBase>();
  params.addClassDescription("Free outflow BC based boundary flux user object for the advection equation using a cell-centered finite volume method.");
  params.addRequiredCoupledVar("u", "Name of the variable to use");
  return params;
}

AEFVFreeOutflowBoundaryFlux::AEFVFreeOutflowBoundaryFlux(const InputParameters & parameters) :
    BoundaryFluxBase(parameters),
    _uc1(coupledValue("u")),
    _u1(getMaterialProperty<Real>("u"))
{
}

AEFVFreeOutflowBoundaryFlux::~AEFVFreeOutflowBoundaryFlux()
{
}

void
AEFVFreeOutflowBoundaryFlux::computeFlux()
{
  unsigned int _qp = 0;

  // assemble the input vectors, which are
  //   the reconstructed linear monomial
  //   extrapolated at side center from the current and neighbor elements
  std::vector<Real> uvec1 = { _u1[_qp] };

  calcFlux(_current_side, _current_elem->id(), uvec1, _normals[_qp], _flux);
}

void
AEFVFreeOutflowBoundaryFlux::calcFlux(unsigned int /*iside*/,
                                      dof_id_type /*ielem*/,
                                      const std::vector<Real> & uvec1,
                                      const RealVectorValue & dwave,
                                      std::vector<Real> & flux)
{
  mooseAssert(uvec1.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");

  // assume the velocity vector is constant, e.g. = (1., 1., 1.)
  RealVectorValue uadv1(1.0, 1.0, 1.0);

  // assign the size of flux vector, e.g. = 1 for the advection equation
  flux.resize(1);

  // finally calculate the flux
  flux[0] = (uadv1 * dwave) * uvec1[0];
}

void
AEFVFreeOutflowBoundaryFlux::computeJacobian()
{
  unsigned int _qp = 0;

  // assemble the input vectors, which are
  //   the constant monomial from the current element
  std::vector<Real> uvec1 = { _uc1[_qp] };

  // calculate the flux
  calcJacobian(_current_side, _current_elem->id(), uvec1, _normals[_qp], _jac1);
}

void
AEFVFreeOutflowBoundaryFlux::calcJacobian(unsigned int /*iside*/,
                                          dof_id_type /*ielem*/,
                                          const std::vector<Real> & libmesh_dbg_var(uvec1),
                                          const RealVectorValue & /*dwave*/,
                                          DenseMatrix<Real> & /*jac1*/)
{
  mooseAssert(uvec1.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");
}
