/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGSystem.h"
#include "FEProblem.h"
#include "TimeIntegrator.h"
#include "RDGAssembleThread.h"
#include "RDGSlopeReconstructionThread.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"

RDGSystem::RDGAssembly::RDGAssembly(RDGSystem & rdg_system)
  : _rdg_system(rdg_system)
{
}

void
RDGSystem::RDGAssembly::assemble()
{
  _rdg_system.assemble();
}


RDGSystem::RDGSystem(FEProblemBase & fe_problem, const std::string & name)
  : NonlinearSystemBase(fe_problem, fe_problem.es().add_system<TransientLinearImplicitSystem>(name), name),
    _rdg_assembly(*this),
    _sys(fe_problem.es().get_system<TransientLinearImplicitSystem>(name)),
    _need_matrix(true)
{
  _sys.attach_assemble_object(_rdg_assembly);
}

RDGSystem::~RDGSystem()
{
}

void
RDGSystem::solve()
{
  // Clear the iteration counters
  _current_l_its.clear();
  _current_nl_its = 0;

  // Initialize the solution vector using a predictor and known values from nodal bcs
  setInitialSolution();

  _time_integrator->solve();
  _time_integrator->postSolve();

  // store info about the solve
  _n_iters = 0;
  _n_linear_iters = 1;

  // TODO: compute final residual
  _final_residual = 0;
}


void
RDGSystem::stopSolve()
{
}

void
RDGSystem::setupFiniteDifferencedPreconditioner()
{
  mooseError2("Finite differencing not available when using rDG.");
}

bool
RDGSystem::converged()
{
  // _console << "converged = " << _sys.linear_solver->get_converged_reason() << std::endl;

  return _sys.linear_solver->get_converged_reason() > 0;
}

void
RDGSystem::assemble()
{
  residualVector(Moose::KT_TIME).zero();
  residualVector(Moose::KT_NONTIME).zero();

  // residual contributions from the domain
  PARALLEL_TRY
  {
    Moose::perf_log.push("RDGSlopeReconstructionThread()", "Execution");

    ConstElemRange & elem_range = *_mesh.getActiveLocalElementRange();
    RDGSlopeReconstructionThread sr(_fe_problem,
                                    _reconstruction_objects,
                                    _limiting_objects);
    Threads::parallel_reduce(elem_range, sr);
    Moose::perf_log.pop("RDGSlopeReconstructionThread()", "Execution");

    Moose::perf_log.push("RDGAssembleThread()", "Execution");
    RDGAssembleThread as(_fe_problem,
                         _need_matrix ? sys().matrix : NULL,
                         _boundary_flux_objects,
                         _internal_side_flux_objects);
    Threads::parallel_reduce(elem_range, as);
    Moose::perf_log.pop("RDGAssembleThread()", "Execution");
  }
  PARALLEL_CATCH;

  sys().matrix->close();

  // gather all contributions to rhs
  residualVector(Moose::KT_TIME).close();
  residualVector(Moose::KT_NONTIME).close();

  sys().rhs->zero();
  _time_integrator->postStep(*sys().rhs);
  sys().rhs->close();

  // FIXME: need to be false, so we form the matrix only once. This is here until we fix libMesh
  _need_matrix = true;
}
