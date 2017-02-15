/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGSystem.h"
#include "FEProblem.h"
#include "TimeIntegrator.h"
#include "AssembleMassMatrixThread.h"
#include "RDGAssembleThread.h"
#include "RDGSlopeReconstructionThread.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"
#include "libmesh/quadrature.h"

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
    _dt(_fe_problem.dt()),
    _need_matrix(true)
{
  _sys.attach_assemble_object(_rdg_assembly);
  _sys.zero_out_matrix_and_rhs = false;
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
RDGSystem::assembleMassMatrix()
{
  std::cerr << "RDGSystem::assembleMassMatrix()" << std::endl;

  const MeshBase & mesh = this->mesh().getMesh();
  const unsigned int dim = mesh.mesh_dimension();

  const DofMap & dof_map = _sys.get_dof_map();
  FEType fe_type(CONSTANT, MONOMIAL);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  UniquePtr<QBase> qrule(QBase::build(QGAUSS, dim, FIRST));
  fe->attach_quadrature_rule(qrule.get());

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  DenseMatrix<Number> Ke;
  Ke.resize(1, 1);
  std::vector<dof_id_type> dof_indices;

  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);

    Ke(0, 0) = JxW[0] * (phi[0][0] * phi[0][0]) / _dt;

    _sys.matrix->add_matrix (Ke, dof_indices);
  }
}

void
RDGSystem::slopeReconstruction()
{
  std::cerr << "RDGSystem::slopeReconstruction()" << std::endl;

  const MeshBase & mesh = this->mesh().getMesh();
  const unsigned int dim = mesh.mesh_dimension();

  const DofMap & dof_map = _sys.get_dof_map();
  FEType fe_type(CONSTANT, MONOMIAL);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  UniquePtr<QBase> qrule(QBase::build(QGAUSS, dim, FIRST));
  fe->attach_quadrature_rule(qrule.get());

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  Densevector<Number> rhhs;
  rhs.resize(1);

  std::vector<dof_id_type> dof_indices;

  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    // fe->reinit(elem);
    //
    // Ke(0, 0) = JxW[0] * (phi[0][0] * phi[0][0]) / _dt;
    //
    // _sys.matrix->add_matrix (Ke, dof_indices);
  }
}

void
RDGSystem::assembleRHS()
{
  const MeshBase & mesh = this->mesh().getMesh();
  const unsigned int dim = mesh.mesh_dimension();

  const DofMap & dof_map = _sys.get_dof_map();
  FEType fe_type(CONSTANT, MONOMIAL);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  UniquePtr<QBase> qrule(QBase::build(QGAUSS, dim, FIRST));
  fe->attach_quadrature_rule(qrule.get());

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  Densevector<Number> rhhs;
  rhs.resize(1);

  std::vector<dof_id_type> dof_indices;

  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    // fe->reinit(elem);
    //
    // Ke(0, 0) = JxW[0] * (phi[0][0] * phi[0][0]) / _dt;
    //
    // _sys.matrix->add_matrix (Ke, dof_indices);
  }
}

void
RDGSystem::assemble()
{
  // mass matrix

  if (_need_matrix)
  {
    PARALLEL_TRY
    {
      Moose::perf_log.push("Assemble mass matrix", "Execution");
      sys().matrix->zero();

      assembleMassMatrix();
      _need_matrix = false;

      sys().matrix->close();
      // sys().matrix->print(std::cerr);
      Moose::perf_log.pop("Assemble mass matrix", "Execution");
    }
    PARALLEL_CATCH;
  }

  // slope reconstruction
  {
    Moose::perf_log.push("Slope reconstruction", "Execution");
    slopeReconstruction();
    Moose::perf_log.pop("Slope reconstruction", "Execution");
  }

  // rDG assemble
  {
    Moose::perf_log.push("rDG assemble", "Execution");
    assembleRHS();
    Moose::perf_log.pop("rDG assemble", "Execution");
  }

  // right hand side
/*
  residualVector(Moose::KT_TIME).zero();
  residualVector(Moose::KT_NONTIME).zero();

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
                         _boundary_flux_objects,
                         _internal_side_flux_objects);
    Threads::parallel_reduce(elem_range, as);
    Moose::perf_log.pop("RDGAssembleThread()", "Execution");
  }
  PARALLEL_CATCH;

  // gather all contributions to rhs
  residualVector(Moose::KT_TIME).close();
  residualVector(Moose::KT_NONTIME).close();

  sys().rhs->zero();
  _time_integrator->postStep(*sys().rhs);
  sys().rhs->close();
*/

  sys().rhs->zero();
  sys().rhs->close();
}
