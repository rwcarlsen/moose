//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEUSEROBJECTSTHREAD_H
#define COMPUTEUSEROBJECTSTHREAD_H

// MOOSE includes
#include "ThreadedElementLoop.h"

#include "libmesh/elem_range.h"

// libMesh forward declarations
namespace libMesh
{
template <typename T>
class NumericVector;
}

/**
 * Class for threaded computation of UserObjects.
 */
class ComputeUserObjectsThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeUserObjectsThread(
      FEProblemBase & problem,
      SystemBase & sys,
      const MooseObjectWarehouse<ElementUserObject> & elemental_user_objects,
      const MooseObjectWarehouse<SideUserObject> & side_user_objects,
      const MooseObjectWarehouse<InternalSideUserObject> & internal_side_user_objects);
  // Splitting Constructor
  ComputeUserObjectsThread(ComputeUserObjectsThread & x, Threads::split);

  virtual ~ComputeUserObjectsThread();

  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void post() override;
  virtual void subdomainChanged() override;

  void join(const ComputeUserObjectsThread & /*y*/);

protected:
  const NumericVector<Number> & _soln;
  _w

      ///@{
      /// Storage for UserObjects (see FEProblemBase::computeUserObjects)
      const MooseObjectWarehouse<ElementUserObject> & _elemental_user_objects;
  const MooseObjectWarehouse<SideUserObject> & _side_user_objects;
  const MooseObjectWarehouse<InternalSideUserObject> & _internal_side_user_objects;
  ///@}
private:
  template <typename T>
  TheWarehouse::Builder querySubdomain(Interfaces iface, std::vector<T> & results)
  {
    _fe_problem.theWarehouse()
        .build()
        .thread(_tid)
        .subdomain(_subdomain)
        .interfaces(iface)
        .queryInto(results);
  }
  template <typename T>
  TheWarehouse::Builder queryBoundary(Interfaces iface, BoundryID bnd, std::vector<T> & results)
  {
    _fe_problem.theWarehouse().build().thread(_tid).boundary(bnd).interfaces(iface).queryInto(
        results);
  }
};

#endif // COMPUTEUSEROBJECTSTHREAD_H
