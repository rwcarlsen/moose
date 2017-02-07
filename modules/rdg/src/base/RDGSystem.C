/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGSystem.h"
#include "FEProblem.h"
#include "TimeIntegrator.h"

RDGSystem::RDGAssembly::RDGAssembly(RDGSystem & rdg_system)
  : _rdg_system(rdg_system)
{
}

void
RDGSystem::RDGAssembly::assemble()
{
  std::cerr << "RDGSystem::RDGAssembly::assemble()" << std::endl;

  TransientExplicitSystem & sys = _rdg_system.sys();
}


RDGSystem::RDGSystem(FEProblemBase & fe_problem, const std::string & name)
  : NonlinearSystemBase(fe_problem, fe_problem.es().add_system<TransientExplicitSystem>(name), name),
    _rdg_assembly(*this),
    _sys(fe_problem.es().get_system<TransientExplicitSystem>(name)),
    _mass_matrix(_sys.add_vector("mass_matrix"))
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

  std::cerr << "_time_integrator = " << _time_integrator << std::endl;

  _time_integrator->solve();
  _time_integrator->postSolve();

  // store info about the solve
  _n_iters = 0;
  // TODO: compute final residual
  _final_residual = 0;

  // FIXME
// #ifdef LIBMESH_HAVE_PETSC
//   _n_linear_iters = static_cast<PetscNonlinearSolver<Real> &>(*_transient_sys.nonlinear_solver).get_total_linear_iterations();
// #endif
}


void
RDGSystem::stopSolve()
{
}

void
RDGSystem::setupFiniteDifferencedPreconditioner()
{
  mooseError("Finite differencing not available when using rDG.");
}

bool
RDGSystem::converged()
{
  return true;
  // return _transient_sys.nonlinear_solver->converged;
}
