/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGSLOPERECONSTRUCTIONTHREAD_H
#define RDGSLOPERECONSTRUCTIONTHREAD_H

#include "ThreadedElementLoop.h"

// libMesh includes
#include "libmesh/elem_range.h"

// Forward declarations
class FEProblemBase;
class NonlinearSystemBase;
class SlopeReconstructionBase;
class SlopeLimitingBase;

class RDGSlopeReconstructionThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  RDGSlopeReconstructionThread(FEProblemBase & fe_problem,
                               const MooseObjectWarehouse<SlopeReconstructionBase> & sr_objects,
                               const MooseObjectWarehouse<SlopeLimitingBase> & limiting_objects);
  // Splitting Constructor
  RDGSlopeReconstructionThread(RDGSlopeReconstructionThread & x, Threads::split split);

  virtual ~RDGSlopeReconstructionThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInterface(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const RDGSlopeReconstructionThread & /*y*/);

protected:
  NonlinearSystemBase & _nl;

  const MooseObjectWarehouse<SlopeReconstructionBase> & _sr_objects;
  const MooseObjectWarehouse<SlopeLimitingBase> & _limiting_objects;
};

#endif //RDGSLOPERECONSTRUCTIONTHREAD_H
