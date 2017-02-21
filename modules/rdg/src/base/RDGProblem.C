/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGProblem.h"
#include "RDGSystem.h"

// because of the addUserObject method
#include "DisplacedProblem.h"
#include "ElementUserObject.h"
#include "NodalUserObject.h"
#include "SideUserObject.h"
#include "InternalSideUserObject.h"
#include "GeneralUserObject.h"

template<>
InputParameters validParams<RDGProblem>()
{
  InputParameters params = validParams<FEProblemBase>();
  return params;
}

RDGProblem::RDGProblem(const InputParameters & parameters)
  : FEProblemBase(parameters)
{
  _nl = _rdg_sys = new RDGSystem(*this, "rdg0");
  _aux = new AuxiliarySystem(*this, "aux0");

  newAssemblyArray(*_nl);

  initNullSpaceVectors(parameters, *_nl);
}

RDGProblem::~RDGProblem()
{
  FEProblemBase::deleteAssemblyArray();

  delete _nl;
  delete _aux;
}

void
RDGProblem::solve()
{
  Moose::perf_log.push("solve()", "Execution");

  _rdg_sys->needMatrix(!_has_jacobian || !_const_jacobian);

  // _nl->computeTimeDerivatives();

  if (_solve)
  {
    _nl->solve();
    _nl->update();
    _has_jacobian = true;
  }

  Moose::perf_log.pop("solve()", "Execution");
}

void
RDGProblem::addUserObject(std::string user_object_name, const std::string & name, InputParameters parameters)
{
  setInputParametersFEProblem(parameters);
  if (_displaced_problem != NULL && parameters.get<bool>("use_displaced_mesh"))
    parameters.set<SubProblem *>("_subproblem") = _displaced_problem.get();
  else
  {
    if (_displaced_problem == NULL && parameters.get<bool>("use_displaced_mesh"))
    {
      // We allow UserObjects to request that they use_displaced_mesh,
      // but then be overridden when no displacements variables are
      // provided in the Mesh block.  If that happened, update the value
      // of use_displaced_mesh appropriately for this UserObject.
      if (parameters.have_parameter<bool>("use_displaced_mesh"))
        parameters.set<bool>("use_displaced_mesh") = false;
    }

    parameters.set<SubProblem *>("_subproblem") = this;
  }

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    // Create the UserObject
    MooseSharedPointer<UserObject> user_object = _factory.create<UserObject>(user_object_name, name, parameters, tid);
    _all_user_objects.addObject(user_object, tid);

    MooseSharedPointer<SlopeReconstructionBase> sruo = MooseSharedNamespace::dynamic_pointer_cast<SlopeReconstructionBase>(user_object);
    MooseSharedPointer<SlopeLimitingBase> sluo = MooseSharedNamespace::dynamic_pointer_cast<SlopeLimitingBase>(user_object);
    MooseSharedPointer<BoundaryFluxBase> bfuo = MooseSharedNamespace::dynamic_pointer_cast<BoundaryFluxBase>(user_object);
    MooseSharedPointer<InternalSideFluxBase> isfuo = MooseSharedNamespace::dynamic_pointer_cast<InternalSideFluxBase>(user_object);

    if (sruo)
      _rdg_sys->_reconstruction_objects.addObject(sruo, tid);
    else if (sluo)
      _rdg_sys->_limiting_objects.addObject(sluo, tid);
    else if (bfuo)
      _rdg_sys->_boundary_flux_objects.addObject(bfuo, tid);
    else if (isfuo)
      _rdg_sys->_internal_side_flux_objects.addObject(isfuo, tid);
    else
    {
      // Attempt to create all the possible UserObject types
      MooseSharedPointer<ElementUserObject> euo = MooseSharedNamespace::dynamic_pointer_cast<ElementUserObject>(user_object);
      MooseSharedPointer<SideUserObject> suo = MooseSharedNamespace::dynamic_pointer_cast<SideUserObject>(user_object);
      MooseSharedPointer<InternalSideUserObject> isuo = MooseSharedNamespace::dynamic_pointer_cast<InternalSideUserObject>(user_object);
      MooseSharedPointer<NodalUserObject> nuo = MooseSharedNamespace::dynamic_pointer_cast<NodalUserObject>(user_object);
      MooseSharedPointer<GeneralUserObject> guo = MooseSharedNamespace::dynamic_pointer_cast<GeneralUserObject>(user_object);

      // Account for displaced mesh use
      if (_displaced_problem != NULL && parameters.get<bool>("use_displaced_mesh"))
      {
        if (euo || nuo)
          _reinit_displaced_elem = true;
        else if (suo)
          _reinit_displaced_face = true;
      }

      // Add the object to the correct warehouse
      if (guo)
      {
        _general_user_objects.addObject(guo);
        break; // not threaded
      }
      else if (nuo)
        _nodal_user_objects.addObject(nuo, tid);
      else if (suo)
        _side_user_objects.addObject(suo, tid);
      else if (isuo)
        _internal_side_user_objects.addObject(isuo, tid);
      else if (euo)
        _elemental_user_objects.addObject(euo, tid);
    }
  }

/*
  FEProblemBase::addUserObject(user_object_name, name, parameters);

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    if (_elemental_user_objects.hasActiveObject(name, tid))
    {
      MooseSharedPointer<UserObject> uo = _elemental_user_objects.getActiveObject(name, tid);

      // put slope reconstruction object into our warehouse
      MooseSharedPointer<SlopeReconstructionBase> sruo = MooseSharedNamespace::dynamic_pointer_cast<SlopeReconstructionBase>(uo);
      if (sruo)
        _rdg_sys->_reconstruction_objects.addObject(sruo, tid);

      // put slope limiting object into our warehouse
      MooseSharedPointer<SlopeLimitingBase> sluo = MooseSharedNamespace::dynamic_pointer_cast<SlopeLimitingBase>(uo);
      if (sluo)
        _rdg_sys->_limiting_objects.addObject(sluo, tid);
    }
*/
}
