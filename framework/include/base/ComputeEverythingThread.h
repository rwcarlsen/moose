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

#ifndef COMPUTEEVERYTHINGTHREAD_H
#define COMPUTEEVERYTHINGTHREAD_H

// libMesh includes
#include "libmesh/elem_range.h"

// MOOSE includes
#include "ThreadedElementLoop.h"
#include "MooseObjectWarehouse.h"
#include "AuxKernel.h"

// Forward declarations
class FEProblemBase;
class AuxiliarySystem;
class ElementUserObject;
class ElementUserObject;
class SystemBase;

namespace libMesh
{
template <typename T>
class NumericVector;
}

class ComputeEverythingThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeEverythingThread(
      FEProblemBase & problem,
      SystemBase & sys,
      const MooseObjectWarehouse<AuxKernel> & storage,
      const MooseObjectWarehouse<ElementUserObject> & pre_elemental_user_objects,
      const MooseObjectWarehouse<ElementUserObject> & post_elemental_user_objects,
      const MooseObjectWarehouse<SideUserObject> & side_user_objects,
      const MooseObjectWarehouse<InternalSideUserObject> & internal_side_user_objects);
  // Splitting Constructor
  ComputeEverythingThread(ComputeEverythingThread & x, Threads::split split);

  virtual ~ComputeEverythingThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void postElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;

  virtual void post() override;

  void join(const ComputeEverythingThread & /*y*/);

protected:
  AuxiliarySystem & _aux_sys;

  /// Storage object containing active AuxKernel objects
  const MooseObjectWarehouse<AuxKernel> & _aux_kernels;

  const NumericVector<Number> & _soln;
  const MooseObjectWarehouse<ElementUserObject> & _pre_elemental_user_objects;
  const MooseObjectWarehouse<ElementUserObject> & _post_elemental_user_objects;
  const MooseObjectWarehouse<SideUserObject> & _side_user_objects;
  const MooseObjectWarehouse<InternalSideUserObject> & _internal_side_user_objects;
};

#endif // COMPUTEELEMAUXVARSTHREAD_H
