/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGAssembleThread.h"
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

RDGAssembleThread::RDGAssembleThread(FEProblemBase & fe_problem,
                                     const MooseObjectWarehouse<BoundaryFluxBase> & bf_objects,
                                     const MooseObjectWarehouse<InternalSideFluxBase> & isf_objects)
  : ThreadedElementLoop<ConstElemRange>(fe_problem),
    _nl(fe_problem.getNonlinearSystemBase()),
    _num_cached(0),
    _integrated_bcs(_nl.getIntegratedBCWarehouse()),
    _dg_kernels(_nl.getDGKernelWarehouse()),
    _interface_kernels(_nl.getInterfaceKernelWarehouse()),
    _kernels(_nl.getKernelWarehouse()),
    _boundary_flux_objects(bf_objects),
    _internal_side_flux_objects(isf_objects)
{
}

// Splitting Constructor
RDGAssembleThread::RDGAssembleThread(RDGAssembleThread & x, Threads::split split) :
    ThreadedElementLoop<ConstElemRange>(x, split),
    _nl(x._nl),
    _num_cached(0),
    _integrated_bcs(x._integrated_bcs),
    _dg_kernels(x._dg_kernels),
    _interface_kernels(x._interface_kernels),
    _kernels(x._kernels),
    _boundary_flux_objects(x._boundary_flux_objects),
    _internal_side_flux_objects(x._internal_side_flux_objects)
{
}

RDGAssembleThread::~RDGAssembleThread()
{
}

void
RDGAssembleThread::subdomainChanged()
{
  _fe_problem.subdomainSetup(_subdomain, _tid);

  // Update variable Dependencies
  std::set<MooseVariable *> needed_moose_vars;
  _kernels.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);
  _integrated_bcs.updateBoundaryVariableDependency(needed_moose_vars, _tid);
  _dg_kernels.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);
  _interface_kernels.updateBoundaryVariableDependency(needed_moose_vars, _tid);

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}

void
RDGAssembleThread::onElement(const Elem * elem)
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
      kernel->computeResidual();
  }

  // const MooseObjectWarehouse<KernelBase> & ntk_warehouse = _nl.getNonTimeKernelWarehouse();
  // if (ntk_warehouse.hasActiveBlockObjects(_subdomain, _tid))
  // {
  //   const std::vector<MooseSharedPointer<KernelBase> > & kernels = ntk_warehouse.getActiveBlockObjects(_subdomain, _tid);
  //   for (const auto & kernel : kernels)
  //     kernel->computeResidual();
  // }
}

void
RDGAssembleThread::onBoundary(const Elem *elem, unsigned int side, BoundaryID bnd_id)
{
  if (_integrated_bcs.hasActiveBoundaryObjects(bnd_id, _tid))
  {
    const std::vector<MooseSharedPointer<IntegratedBC> > & bcs = _integrated_bcs.getActiveBoundaryObjects(bnd_id, _tid);

    _fe_problem.reinitElemFace(elem, side, bnd_id, _tid);

    // Set up Sentinel class so that, even if reinitMaterialsFace() throws, we
    // still remember to swap back during stack unwinding.
    SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);

    _fe_problem.reinitMaterialsFace(elem->subdomain_id(), _tid);
    _fe_problem.reinitMaterialsBoundary(bnd_id, _tid);

    // Set the active boundary id so that BoundaryRestrictable::_boundary_id is correct
    _fe_problem.setCurrentBoundaryID(bnd_id);

    // boundary fluxes
    if (_boundary_flux_objects.hasActiveBoundaryObjects(bnd_id, _tid))
    {
      const std::vector<MooseSharedPointer<BoundaryFluxBase> > & objects = _boundary_flux_objects.getActiveBoundaryObjects(_subdomain, _tid);
      for (const auto & uo : objects)
        uo->computeFlux();

      // TODO: compute Jacobian if needed
    }

    for (const auto & bc : bcs)
    {
      if (bc->shouldApply())
        bc->computeResidual();
    }

    // Set active boundary id to invalid
    _fe_problem.setCurrentBoundaryID(Moose::INVALID_BOUNDARY_ID);
  }
}

void
RDGAssembleThread::onInterface(const Elem * /*elem*/, unsigned int /*side*/, BoundaryID bnd_id)
{
  if (_interface_kernels.hasActiveBoundaryObjects(bnd_id, _tid))
  {
#if 0
    // Pointer to the neighbor we are currently working on.
    const Elem * neighbor = elem->neighbor(side);

    if (!(neighbor->level() == elem->level()))
      mooseError2("Sorry, interface kernels do not work with mesh adaptivity");

    if (neighbor->active())
    {
      _fe_problem.reinitNeighbor(elem, side, _tid);

      // Set up Sentinels so that, even if one of the reinitMaterialsXXX() calls throws, we
      // still remember to swap back during stack unwinding.
      SwapBackSentinel face_sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);
      _fe_problem.reinitMaterialsFace(elem->subdomain_id(), _tid);

      SwapBackSentinel neighbor_sentinel(_fe_problem, &FEProblem::swapBackMaterialsNeighbor, _tid);
      _fe_problem.reinitMaterialsNeighbor(neighbor->subdomain_id(), _tid);

      const std::vector<MooseSharedPointer<InterfaceKernel> > & int_ks = _interface_kernels.getActiveBoundaryObjects(bnd_id, _tid);
      for (const auto & interface_kernel : int_ks)
        interface_kernel->computeResidual();

      {
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
        _fe_problem.addResidualNeighbor(_tid);
      }
    }
#endif
  }
}

void
RDGAssembleThread::onInternalSide(const Elem * elem, unsigned int side)
{
#if 1
  if (_dg_kernels.hasActiveBlockObjects(_subdomain, _tid))
  {
    // Pointer to the neighbor we are currently working on.
    const Elem * neighbor = elem->neighbor(side);

    // Get the global id of the element and the neighbor
    const dof_id_type
      elem_id = elem->id(),
      neighbor_id = neighbor->id();

    if ((neighbor->active() && (neighbor->level() == elem->level()) && (elem_id < neighbor_id)) || (neighbor->level() < elem->level()))
    {
      _fe_problem.reinitNeighbor(elem, side, _tid);

      // Set up Sentinels so that, even if one of the reinitMaterialsXXX() calls throws, we
      // still remember to swap back during stack unwinding.
      SwapBackSentinel face_sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);
      _fe_problem.reinitMaterialsFace(elem->subdomain_id(), _tid);

      SwapBackSentinel neighbor_sentinel(_fe_problem, &FEProblem::swapBackMaterialsNeighbor, _tid);
      _fe_problem.reinitMaterialsNeighbor(neighbor->subdomain_id(), _tid);

      // internal side fluxes
      if (_internal_side_flux_objects.hasActiveBlockObjects(_subdomain, _tid))
      {
        const std::vector<MooseSharedPointer<InternalSideFluxBase> > & objects = _internal_side_flux_objects.getActiveBlockObjects(_subdomain, _tid);
        for (const auto & uo : objects)
          uo->computeFlux();

        // TODO: compute Jacobian if needed
      }

      const std::vector<MooseSharedPointer<DGKernel> > & dgks = _dg_kernels.getActiveBlockObjects(_subdomain, _tid);
      for (const auto & dg_kernel : dgks)
        if (dg_kernel->hasBlocks(neighbor->subdomain_id()))
          dg_kernel->computeResidual();

      {
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
        _fe_problem.addResidualNeighbor(_tid);
      }
    }
  }
#endif
}

void
RDGAssembleThread::postElement(const Elem * /*elem*/)
{
  _fe_problem.cacheResidual(_tid);
  _num_cached++;

  if (_num_cached % 20 == 0)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    _fe_problem.addCachedResidual(_tid);
  }
}

void
RDGAssembleThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
  _fe_problem.addCachedResidual(_tid);
}


void
RDGAssembleThread::join(const RDGAssembleThread & /*y*/)
{
}
