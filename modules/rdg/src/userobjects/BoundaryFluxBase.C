/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BoundaryFluxBase.h"

// Static mutex definition
Threads::spin_mutex BoundaryFluxBase::_mutex;

template<>
InputParameters validParams<BoundaryFluxBase>()
{
  InputParameters params = validParams<SideUserObject>();
  return params;
}

BoundaryFluxBase::BoundaryFluxBase(const InputParameters & parameters) :
    SideUserObject(parameters)
{
}

void
BoundaryFluxBase::initialize()
{
}

void
BoundaryFluxBase::execute()
{
}

void
BoundaryFluxBase::finalize()
{
}

void
BoundaryFluxBase::threadJoin(const UserObject & uo)
{
}

void
BoundaryFluxBase::computeFlux()
{
}

void
BoundaryFluxBase::computeJacobian()
{
}

const std::vector<Real> &
BoundaryFluxBase::getFlux() const
{
  return _flux;
}

const DenseMatrix<Real> &
BoundaryFluxBase::getJacobian() const
{
  return _jac1;
}

/*
const std::vector<Real> &
BoundaryFluxBase::getFlux(unsigned int iside,
                          dof_id_type ielem,
                          const std::vector<Real> & uvec1,
                          const RealVectorValue & dwave,
                          THREAD_ID tid) const
{
  Threads::spin_mutex::scoped_lock lock(_mutex);
  if (_cached_elem_id != ielem || _cached_side_id != iside)
  {
    _cached_elem_id = ielem;
    _cached_side_id = iside;

    calcFlux(iside,
             ielem,
             uvec1,
             dwave,
             _flux[tid]);
  }
  return _flux[tid];
}

const DenseMatrix<Real> &
BoundaryFluxBase::getJacobian(unsigned int iside,
                              dof_id_type ielem,
                              const std::vector<Real> & uvec1,
                              const RealVectorValue & dwave,
                              THREAD_ID tid) const
{
  Threads::spin_mutex::scoped_lock lock(_mutex);
  if (_cached_elem_id != ielem || _cached_side_id != iside)
  {
    _cached_elem_id = ielem;
    _cached_side_id = iside;

    calcJacobian(iside,
                 ielem,
                 uvec1,
                 dwave,
                 _jac1[tid]);
  }
  return _jac1[tid];
}
*/
