/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "InternalSideFluxBase.h"

// Static mutex definition
Threads::spin_mutex InternalSideFluxBase::_mutex;

template<>
InputParameters validParams<InternalSideFluxBase>()
{
  InputParameters params = validParams<InternalSideUserObject>();
  params.addClassDescription("A base class for computing and caching internal side flux.");
  return params;
}

InternalSideFluxBase::InternalSideFluxBase(const InputParameters & parameters) :
    InternalSideUserObject(parameters)
{
}

void
InternalSideFluxBase::initialize()
{
}

void
InternalSideFluxBase::execute()
{
}

void
InternalSideFluxBase::finalize()
{
}

void
InternalSideFluxBase::threadJoin(const UserObject & uo)
{
}

void
InternalSideFluxBase::computeFlux()
{
}

void
InternalSideFluxBase::computeJacobian()
{
}

const std::vector<Real> &
InternalSideFluxBase::getFlux() const
{
  return _flux;
}

const DenseMatrix<Real> &
InternalSideFluxBase::getJacobian(Moose::DGResidualType type) const
{
  if (type == Moose::Element)
    return _jac1;
  else
    return _jac2;
}
