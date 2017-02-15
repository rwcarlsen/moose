/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef ASSEMBLEMASSMATRIXTHREAD_H
#define ASSEMBLEMASSMATRIXTHREAD_H

#include "ThreadedElementLoop.h"

// libMesh includes
#include "libmesh/elem_range.h"

// Forward declarations
class FEProblemBase;
class NonlinearSystemBase;
class KernelWarehouse;

class AssembleMassMatrixThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  AssembleMassMatrixThread(FEProblemBase & fe_problem, SparseMatrix<Number> & matrix);
  AssembleMassMatrixThread(AssembleMassMatrixThread & x, Threads::split split);

  virtual ~AssembleMassMatrixThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInterface(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const AssembleMassMatrixThread & /*y*/);

protected:
  NonlinearSystemBase & _nl;
  SparseMatrix<Number> & _matrix;
  unsigned int _num_cached;
  const KernelWarehouse & _kernels;
};

#endif //ASSEMBLEMASSMATRIXTHREAD_H
