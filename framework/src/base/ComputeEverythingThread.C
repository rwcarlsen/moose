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

#include "ComputeEverythingThread.h"

// MOOSE includes
#include "AuxiliarySystem.h"
#include "AuxKernel.h"
#include "SwapBackSentinel.h"
#include "FEProblem.h"
#include "SystemBase.h"
#include "ElementUserObject.h"
#include "ShapeElementUserObject.h"
#include "Problem.h"

// libMesh includes
#include "libmesh/threads.h"
#include "libmesh/numeric_vector.h"

ComputeEverythingThread::ComputeEverythingThread(
      FEProblemBase & problem,
      SystemBase & sys,
      const MooseObjectWarehouse<AuxKernel> & storage,
      const MooseObjectWarehouse<ElementUserObject> & pre_elemental_user_objects,
      const MooseObjectWarehouse<ElementUserObject> & post_elemental_user_objects,
      const MooseObjectWarehouse<SideUserObject> & pre_side_user_objects,
      const MooseObjectWarehouse<SideUserObject> & post_side_user_objects,
      const MooseObjectWarehouse<InternalSideUserObject> & pre_internal_side_user_objects,
      const MooseObjectWarehouse<InternalSideUserObject> & post_internal_side_user_objects)
  : ThreadedElementLoop<ConstElemRange>(problem),
    _aux_sys(problem.getAuxiliarySystem()),
    _aux_kernels(storage),
    _soln(*sys.currentSolution()),
    _pre_elemental_user_objects(pre_elemental_user_objects),
    _post_elemental_user_objects(post_elemental_user_objects),
    _pre_side_user_objects(pre_side_user_objects),
    _post_side_user_objects(post_side_user_objects),
    _pre_internal_side_user_objects(pre_internal_side_user_objects),
    _post_internal_side_user_objects(post_internal_side_user_objects)
{
}

// Splitting Constructor
ComputeEverythingThread::ComputeEverythingThread(ComputeEverythingThread & x,
                                                   Threads::split /*split*/)
  : ThreadedElementLoop<ConstElemRange>(x._fe_problem),
    _aux_sys(x._aux_sys),
    _aux_kernels(x._aux_kernels),
    _soln(x._soln),
    _pre_elemental_user_objects(x._pre_elemental_user_objects),
    _post_elemental_user_objects(x._post_elemental_user_objects),
    _pre_side_user_objects(x._pre_side_user_objects),
    _post_side_user_objects(x._post_side_user_objects),
    _pre_internal_side_user_objects(x._pre_internal_side_user_objects),
    _post_internal_side_user_objects(x._post_internal_side_user_objects)
{
}

ComputeEverythingThread::~ComputeEverythingThread() {}

void
ComputeEverythingThread::subdomainChanged()
{
  ////////////// ElemAuxVars /////////////
  for (const auto & it : _aux_sys._elem_vars[_tid])
  {
    MooseVariable * var = it.second;
    var->prepareAux();
  }

  std::set<MooseVariable *> needed_moose_vars;
  std::set<unsigned int> needed_mat_props;

  if (_aux_kernels.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<std::shared_ptr<AuxKernel>> & kernels =
        _aux_kernels.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & aux : kernels)
    {
      aux->subdomainSetup();
      const std::set<MooseVariable *> & mv_deps = aux->getMooseVariableDependencies();
      const std::set<unsigned int> & mp_deps = aux->getMatPropDependencies();
      needed_moose_vars.insert(mv_deps.begin(), mv_deps.end());
      needed_mat_props.insert(mp_deps.begin(), mp_deps.end());
    }
  }

  ////////////// UserObjects /////////////
  _pre_elemental_user_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);
  _pre_elemental_user_objects.updateBlockMatPropDependency(_subdomain, needed_mat_props, _tid);
  _pre_elemental_user_objects.subdomainSetup(_subdomain, _tid);
  _post_elemental_user_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);
  _post_elemental_user_objects.updateBlockMatPropDependency(_subdomain, needed_mat_props, _tid);
  _post_elemental_user_objects.subdomainSetup(_subdomain, _tid);

  _pre_side_user_objects.updateBoundaryVariableDependency(needed_moose_vars, _tid);
  _post_side_user_objects.updateBoundaryVariableDependency(needed_moose_vars, _tid);
  _pre_internal_side_user_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);
  _post_internal_side_user_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);

  _pre_side_user_objects.updateBoundaryMatPropDependency(needed_mat_props, _tid);
  _post_side_user_objects.updateBoundaryMatPropDependency(needed_mat_props, _tid);
  _pre_internal_side_user_objects.updateBlockMatPropDependency(_subdomain, needed_mat_props, _tid);
  _post_internal_side_user_objects.updateBlockMatPropDependency(_subdomain, needed_mat_props, _tid);

  _pre_side_user_objects.subdomainSetup(_tid);
  _post_side_user_objects.subdomainSetup(_tid);
  _pre_internal_side_user_objects.subdomainSetup(_subdomain, _tid);
  _post_internal_side_user_objects.subdomainSetup(_subdomain, _tid);

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.setActiveMaterialProperties(needed_mat_props, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}

void
ComputeEverythingThread::onElement(const Elem * elem)
{
  _fe_problem.prepare(elem, _tid);
  _fe_problem.reinitElem(elem, _tid);

  // Set up the sentinel so that, even if reinitMaterials() throws, we
  // still remember to swap back.
  SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterials, _tid, true);
  _fe_problem.reinitMaterials(elem->subdomain_id(), _tid);

  /////////// UserObjects ///////////////
  if (_pre_elemental_user_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    const auto & objects = _pre_elemental_user_objects.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & uo : objects)
      uo->execute();
  }

  // UserObject Jacobians
  if (_fe_problem.currentlyComputingJacobian() &&
      _pre_elemental_user_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    // Prepare shape functions for ShapeElementUserObjects
    std::vector<MooseVariable *> jacobian_moose_vars =
        _fe_problem.getUserObjectJacobianVariables(_tid);
    for (auto & jvar : jacobian_moose_vars)
    {
      unsigned int jvar_id = jvar->number();
      std::vector<dof_id_type> & dof_indices = jvar->dofIndices();

      _fe_problem.prepareShapes(jvar_id, _tid);

      const auto & e_objects = _pre_elemental_user_objects.getActiveBlockObjects(_subdomain, _tid);
      for (const auto & uo : e_objects)
      {
        auto shape_element_uo = std::dynamic_pointer_cast<ShapeElementUserObject>(uo);
        if (shape_element_uo)
          shape_element_uo->executeJacobianWrapper(jvar_id, dof_indices);
      }
    }
  }

  ////////////// ElemAuxVars /////////////
  if (_aux_kernels.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<std::shared_ptr<AuxKernel>> & kernels =
        _aux_kernels.getActiveBlockObjects(_subdomain, _tid);

    for (const auto & aux : kernels)
      aux->compute();

    // update the solution vector
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & it : _aux_sys._elem_vars[_tid])
      {
        MooseVariable * var = it.second;
        var->insert(_aux_sys.solution());
      }
    }
  }

  /////////// UserObjects ///////////////
  if (_post_elemental_user_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    const auto & objects = _post_elemental_user_objects.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & uo : objects)
      uo->execute();
  }

  // UserObject Jacobians
  if (_fe_problem.currentlyComputingJacobian() &&
      _post_elemental_user_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    // Prepare shape functions for ShapeElementUserObjects
    std::vector<MooseVariable *> jacobian_moose_vars =
        _fe_problem.getUserObjectJacobianVariables(_tid);
    for (auto & jvar : jacobian_moose_vars)
    {
      unsigned int jvar_id = jvar->number();
      std::vector<dof_id_type> & dof_indices = jvar->dofIndices();

      _fe_problem.prepareShapes(jvar_id, _tid);

      const auto & e_objects = _post_elemental_user_objects.getActiveBlockObjects(_subdomain, _tid);
      for (const auto & uo : e_objects)
      {
        auto shape_element_uo = std::dynamic_pointer_cast<ShapeElementUserObject>(uo);
        if (shape_element_uo)
          shape_element_uo->executeJacobianWrapper(jvar_id, dof_indices);
      }
    }
  }
}

void
ComputeEverythingThread::onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id)
{
  _fe_problem.reinitElemFace(elem, side, bnd_id, _tid);

  // Set up Sentinel class so that, even if reinitMaterialsFace() throws, we
  // still remember to swap back during stack unwinding.
  SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);
  _fe_problem.reinitMaterialsFace(_subdomain, _tid);
  _fe_problem.reinitMaterialsBoundary(bnd_id, _tid);

  _fe_problem.setCurrentBoundaryID(bnd_id);

  //////////// PRE_AUX UserObjects
  if (_pre_side_user_objects.hasActiveBoundaryObjects(bnd_id, _tid))
  {
    const auto & objects = _pre_side_user_objects.getActiveBoundaryObjects(bnd_id, _tid);
    for (const auto & uo : objects)
      uo->execute();

    // UserObject Jacobians
    if (_fe_problem.currentlyComputingJacobian())
    {
      // Prepare shape functions for ShapeSideUserObjects
      std::vector<MooseVariable *> jacobian_moose_vars =
          _fe_problem.getUserObjectJacobianVariables(_tid);
      for (auto & jvar : jacobian_moose_vars)
      {
        unsigned int jvar_id = jvar->number();
        std::vector<dof_id_type> & dof_indices = jvar->dofIndices();

        _fe_problem.prepareFaceShapes(jvar_id, _tid);

        for (const auto & uo : objects)
        {
          auto shape_side_uo = std::dynamic_pointer_cast<ShapeSideUserObject>(uo);
          if (shape_side_uo)
            shape_side_uo->executeJacobianWrapper(jvar_id, dof_indices);
        }
      }
    }

    _fe_problem.setCurrentBoundaryID(Moose::INVALID_BOUNDARY_ID);
  }

  //////////// POST_AUX UserObjects
  if (_post_side_user_objects.hasActiveBoundaryObjects(bnd_id, _tid))
  {
    const auto & objects = _post_side_user_objects.getActiveBoundaryObjects(bnd_id, _tid);
    for (const auto & uo : objects)
      uo->execute();

    // UserObject Jacobians
    if (_fe_problem.currentlyComputingJacobian())
    {
      // Prepare shape functions for ShapeSideUserObjects
      std::vector<MooseVariable *> jacobian_moose_vars =
          _fe_problem.getUserObjectJacobianVariables(_tid);
      for (auto & jvar : jacobian_moose_vars)
      {
        unsigned int jvar_id = jvar->number();
        std::vector<dof_id_type> & dof_indices = jvar->dofIndices();

        _fe_problem.prepareFaceShapes(jvar_id, _tid);

        for (const auto & uo : objects)
        {
          auto shape_side_uo = std::dynamic_pointer_cast<ShapeSideUserObject>(uo);
          if (shape_side_uo)
            shape_side_uo->executeJacobianWrapper(jvar_id, dof_indices);
        }
      }
    }

    _fe_problem.setCurrentBoundaryID(Moose::INVALID_BOUNDARY_ID);
  }
}

void
ComputeUserObjectsThread::onInternalSide(const Elem * elem, unsigned int side)
{
  // Pointer to the neighbor we are currently working on.
  const Elem * neighbor = elem->neighbor_ptr(side);

  // Get the global id of the element and the neighbor
  const dof_id_type elem_id = elem->id(), neighbor_id = neighbor->id();

  if (!_internal_side_user_objects.hasActiveBlockObjects(_subdomain, _tid))
    return;
  if (!((neighbor->active() && (neighbor->level() == elem->level()) && (elem_id < neighbor_id)) ||
        (neighbor->level() < elem->level())))
    return;

  _fe_problem.prepareFace(elem, _tid);
  _fe_problem.reinitNeighbor(elem, side, _tid);

  // Set up Sentinels so that, even if one of the reinitMaterialsXXX() calls throws, we
  // still remember to swap back during stack unwinding.
  SwapBackSentinel face_sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);
  _fe_problem.reinitMaterialsFace(elem->subdomain_id(), _tid);

  SwapBackSentinel neighbor_sentinel(_fe_problem, &FEProblem::swapBackMaterialsNeighbor, _tid);
  _fe_problem.reinitMaterialsNeighbor(neighbor->subdomain_id(), _tid);

  const auto & objects = _internal_side_user_objects.getActiveBlockObjects(_subdomain, _tid);
  for (const auto & uo : objects)
  {
    if (!uo->blockRestricted())
      uo->execute();
    else if (uo->hasBlocks(neighbor->subdomain_id()))
      uo->execute();
  }
}


void
ComputeEverythingThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
  _fe_problem.clearActiveMaterialProperties(_tid);
}

void
ComputeEverythingThread::join(const ComputeEverythingThread & /*y*/)
{
}
