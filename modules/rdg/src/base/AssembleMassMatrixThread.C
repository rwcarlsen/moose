/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "AssembleMassMatrixThread.h"
#include "NonlinearSystem.h"
#include "Problem.h"
#include "FEProblem.h"
#include "KernelBase.h"
#include "IntegratedBC.h"
#include "DGKernel.h"
#include "InterfaceKernel.h"
#include "Material.h"
#include "TimeKernel.h"
#include "KernelWarehouse.h"
#include "SwapBackSentinel.h"
#include "BoundaryFluxBase.h"
#include "InternalSideFluxBase.h"

// libmesh includes
#include "libmesh/threads.h"

AssembleMassMatrixThread::AssembleMassMatrixThread(FEProblemBase & fe_problem, SparseMatrix<Number> & matrix)
  : ThreadedElementLoop<ConstElemRange>(fe_problem),
    _nl(fe_problem.getNonlinearSystemBase()),
    _matrix(matrix),
    _num_cached(0),
    _kernels(_nl.getKernelWarehouse())
{
}

// Splitting Constructor
AssembleMassMatrixThread::AssembleMassMatrixThread(AssembleMassMatrixThread & x, Threads::split split) :
    ThreadedElementLoop<ConstElemRange>(x, split),
    _nl(x._nl),
    _matrix(x._matrix),
    _num_cached(0),
    _kernels(x._kernels)
{
}

AssembleMassMatrixThread::~AssembleMassMatrixThread()
{
}

void
AssembleMassMatrixThread::subdomainChanged()
{
  _fe_problem.subdomainSetup(_subdomain, _tid);

  // Update variable Dependencies
  std::set<MooseVariable *> needed_moose_vars;
  _kernels.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}

void
AssembleMassMatrixThread::onElement(const Elem * elem)
{
  _fe_problem.prepare(elem, _tid);
  _fe_problem.reinitElem(elem, _tid);

  // Set up Sentinel class so that, even if reinitMaterials() throws, we
  // still remember to swap back during stack unwinding.
  SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterials, _tid);

  _fe_problem.reinitMaterials(_subdomain, _tid);

  const MooseObjectWarehouse<KernelBase> & tk_warehouse = _nl.getTimeKernelWarehouse();
  if (tk_warehouse.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<MooseSharedPointer<KernelBase> > & kernels = tk_warehouse.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & kernel : kernels)
    {
      kernel->subProblem().prepareShapes(kernel->variable().number(), _tid);
      kernel->computeJacobian();
    }
  }
}

void
AssembleMassMatrixThread::onBoundary(const Elem *elem, unsigned int side, BoundaryID bnd_id)
{
}

void
AssembleMassMatrixThread::onInterface(const Elem * /*elem*/, unsigned int /*side*/, BoundaryID bnd_id)
{
}

void
AssembleMassMatrixThread::onInternalSide(const Elem * elem, unsigned int side)
{
}

void
AssembleMassMatrixThread::postElement(const Elem * /*elem*/)
{
  _fe_problem.cacheJacobian(_tid);
  _num_cached++;

  if (_num_cached % 20 == 0)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    _fe_problem.addCachedJacobian(_matrix, _tid);
  }
}

void
AssembleMassMatrixThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
  _fe_problem.addCachedJacobian(_matrix, _tid);
}


void
AssembleMassMatrixThread::join(const AssembleMassMatrixThread & /*y*/)
{
}
