/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGASSEMBLETHREAD_H
#define RDGASSEMBLETHREAD_H

#include "ThreadedElementLoop.h"

// libMesh includes
#include "libmesh/elem_range.h"

// Forward declarations
class FEProblemBase;
class NonlinearSystemBase;
class IntegratedBC;
class DGKernel;
class InterfaceKernel;
class TimeKernel;
class KernelBase;
class KernelWarehouse;
class BoundaryFluxBase;
class InternalSideFluxBase;

class RDGAssembleThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  RDGAssembleThread(FEProblemBase & fe_problem, SparseMatrix<Number> * matrix,
                    const MooseObjectWarehouse<BoundaryFluxBase> & bf_objects,
                    const MooseObjectWarehouse<InternalSideFluxBase> & isf_objects);
  // Splitting Constructor
  RDGAssembleThread(RDGAssembleThread & x, Threads::split split);

  virtual ~RDGAssembleThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInterface(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const RDGAssembleThread & /*y*/);

protected:
  NonlinearSystemBase & _nl;
  SparseMatrix<Number> * _matrix;
  unsigned int _num_cached;

  const MooseObjectWarehouse<IntegratedBC> & _integrated_bcs;
  const MooseObjectWarehouse<DGKernel> & _dg_kernels;
  const MooseObjectWarehouse<InterfaceKernel> & _interface_kernels;
  const KernelWarehouse & _kernels;

  const MooseObjectWarehouse<BoundaryFluxBase> & _boundary_flux_objects;
  const MooseObjectWarehouse<InternalSideFluxBase> & _internal_side_flux_objects;
};

#endif //RDGASSEMBLETHREAD_H
