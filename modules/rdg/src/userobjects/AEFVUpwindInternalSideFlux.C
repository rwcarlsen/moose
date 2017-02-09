/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "AEFVUpwindInternalSideFlux.h"

template<>
InputParameters validParams<AEFVUpwindInternalSideFlux>()
{
  InputParameters params = validParams<InternalSideFluxBase>();
  params.addClassDescription("Upwind numerical flux scheme for the advection equation using a cell-centered finite volume method.");
  params.addRequiredCoupledVar("u", "Name of the variable to use");
  return params;
}

AEFVUpwindInternalSideFlux::AEFVUpwindInternalSideFlux(const InputParameters & parameters) :
    InternalSideFluxBase(parameters),
    _uc1(coupledValue("u")),
    _uc2(coupledNeighborValue("u")),
    _u1(getMaterialProperty<Real>("u")),
    _u2(getNeighborMaterialProperty<Real>("u"))
{
}

AEFVUpwindInternalSideFlux::~AEFVUpwindInternalSideFlux()
{
}

void
AEFVUpwindInternalSideFlux::computeFlux()
{
  std::cerr << "AEFVUpwindInternalSideFlux::computeFlux()" << std::endl;

  unsigned int _qp = 0;

  // assemble the input vectors, which are
  //   the reconstructed linear monomial
  //   extrapolated at side center from the current and neighbor elements
  std::vector<Real> uvec1 = { _u1[_qp] };
  std::vector<Real> uvec2 = { _u2[_qp] };

  calcFlux(_current_side, _current_elem->id(), _neighbor_elem->id(),
           uvec1, uvec2, _normals[_qp], _flux);
  std::cerr << "flux = " << _flux[0] << std::endl;
}


void
AEFVUpwindInternalSideFlux::calcFlux(unsigned int /*iside*/,
                                     unsigned int /*ielem*/,
                                     unsigned int /*ineig*/,
                                     const std::vector<Real> & uvec1,
                                     const std::vector<Real> & uvec2,
                                     const RealVectorValue & dwave,
                                     std::vector<Real> & flux)
{
  mooseAssert(uvec1.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");
  mooseAssert(uvec2.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");

  // assign the size of flux vector, e.g. = 1 for the advection equation
  flux.resize(1);

  // assume a constant velocity on the left
  RealVectorValue uadv1(1.0, 1.0, 1.0);

  // assume a constant velocity on the right
  RealVectorValue uadv2(1.0, 1.0, 1.0);

  // normal velocity on the left and right
  Real vdon1 = uadv1 * dwave;
  Real vdon2 = uadv2 * dwave;

  // calculate the so-called a^plus and a^minus
  Real aplus = 0.5 * (vdon1 + std::abs(vdon1));
  Real amins = 0.5 * (vdon2 - std::abs(vdon2));

  // finally calculate the flux
  flux[0] = aplus * uvec1[0] + amins * uvec2[0];
}

void
AEFVUpwindInternalSideFlux::computeJacobian()
{
  unsigned int _qp = 0;

  // assemble the input vectors, which are
  //   the constant monomial from the current and neighbor elements
  std::vector<Real> uvec1 = { _uc1[_qp] };
  std::vector<Real> uvec2 = { _uc2[_qp] };

  calcJacobian(_current_side, _current_elem->id(), _neighbor_elem->id(),
               uvec1, uvec2, _normals[_qp], _jac1, _jac2);
}

void
AEFVUpwindInternalSideFlux::calcJacobian(unsigned int /*iside*/,
                                         unsigned int /*ielem*/,
                                         unsigned int /*ineig*/,
                                         const std::vector<Real> & libmesh_dbg_var(uvec1),
                                         const std::vector<Real> & libmesh_dbg_var(uvec2),
                                         const RealVectorValue & dwave,
                                         DenseMatrix<Real> & jac1,
                                         DenseMatrix<Real> & jac2)
{
  mooseAssert(uvec1.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");
  mooseAssert(uvec2.size() == 1, "Invalid size for uvec1. Must be single variable coupling.");

  // assign the size of Jacobian matrix, e.g. = (1, 1) for the advection equation
  jac1.resize(1, 1);
  jac2.resize(1, 1);

  // assume a constant velocity on the left
  RealVectorValue uadv1(1.0, 1.0, 1.0);

  // assume a constant velocity on the right
  RealVectorValue uadv2(1.0, 1.0, 1.0);

  // normal velocity on the left and right
  Real vdon1 = uadv1 * dwave;
  Real vdon2 = uadv2 * dwave;

  // calculate the so-called a^plus and a^minus
  Real aplus = 0.5 * (vdon1 + std::abs(vdon1));
  Real amins = 0.5 * (vdon2 - std::abs(vdon2));

  // finally calculate the Jacobian matrix
  jac1(0, 0) = aplus;
  jac2(0, 0) = amins;
}
