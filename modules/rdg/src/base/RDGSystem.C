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
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
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
    _dof_map(_sys.get_dof_map()),
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

  // std::cerr << "rhs = " << std::endl;
  // _sys.rhs->print(std::cerr);
  //
  // std::cerr << "sln =" << std::endl;
  // solution().print(std::cerr);

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
  // std::cerr << "RDGSystem::assembleMassMatrix()" << std::endl;

  const MeshBase & mesh = this->mesh().getMesh();
  const unsigned int dim = mesh.mesh_dimension();

  FEType fe_type(CONSTANT, MONOMIAL);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  UniquePtr<QBase> qrule(QBase::build(QGAUSS, dim, FIRST));
  fe->attach_quadrature_rule(qrule.get());

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();


  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  UniquePtr<QBase>  qface (QBase::build(QGAUSS, dim - 1, FIRST));
  fe_face->attach_quadrature_rule (qface.get());

  const std::vector<Point> & qface_point = fe_face->get_xyz();
  const std::vector<Point> & normals = fe_face->get_normals();


  DenseMatrix<Number> Ke;
  Ke.resize(1, 1);
  std::vector<dof_id_type> dof_indices;

  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;

    _dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);

    Ke(0, 0) = JxW[0] * (phi[0][0] * phi[0][0]) / _dt;
    // std::cerr << "Ke = " << Ke(0, 0) << std::endl;

    dof_id_type elem_id = elem->id();
    _elem_JxW[elem_id] = JxW[0];
    _elem_centroid[elem_id] = elem->centroid();
    _dof_indices[elem_id] = dof_indices;

    _sys.matrix->add_matrix (Ke, dof_indices);

    for (unsigned int side = 0; side<elem->n_sides(); side++)
    {
      std::vector<BoundaryID> boundary_ids = _mesh.getBoundaryIDs(elem, side);

      fe_face->reinit(elem, side);

      std::pair<dof_id_type, unsigned int> key(elem_id, side);
      _face_point[key] = qface_point[0];
      _normals[key] = normals[0];
    }
  }
}


std::vector<RealGradient>
RDGSystem::limitElementSlope(const Elem * elem)
{
  dof_id_type elem_id = elem->id();

  const std::vector<dof_id_type> & dof_indices = _dof_indices[elem_id];
  // _dof_map.dof_indices(elem, dof_indices);


  // you should know how many equations you are solving and assign this number
  // e.g. = 1 (for the advection equation)
  unsigned int nvars = 1;

  // index of conserved variable
  unsigned int iv;

  // index of sides around an element
  unsigned int is;

  // index of the neighbor element
  unsigned int in;

  // number of sides surrounding an element = 2 in 1D
  unsigned int nside = elem->n_sides();

  // number of reconstruction stencils = 3 in 1D
  unsigned int nsten = nside + 1;

  // vector for the gradients of primitive variables
  std::vector<RealGradient> ugrad(nvars, RealGradient(0., 0., 0.));

  // array to store center coordinates of this cell and its neighbor cells
  std::vector<Real> xc(nsten, 0.);

  // the first always stores the current cell
  xc[0] = _elem_centroid[elem_id](0);

  // array for the cell-average values in the current cell and its neighbors
  std::vector <std::vector<Real> > ucell(nsten, std::vector<Real>(nvars, 0.));

  // central slope:
  //                  u_{i+1} - u {i-1}
  //                  -----------------
  //                  x_{i+1} - x_{i-1}

  // left-side slope:
  //                  u_{i} - u {i-1}
  //                  ---------------
  //                  x_{i} - x_{i-1}

  // right-side slope:
  //                  u_{i+1} - u {i}
  //                  ---------------
  //                  x_{i+1} - x_{i}

  // array to store the central and one-sided slopes, where the first should be central slope
  std::vector <std::vector<Real> > sigma(nsten, std::vector<Real>(nvars, 0.));

  // get the cell-average variable in the central cell
  ucell[0][0] = (*_sys.solution)(dof_indices[0]);
  // std::cerr << "ucell[0][0] = " << ucell[0][0] << std::endl;

  // a flag to indicate the boundary side of the current cell

  unsigned int bflag = 0;

  // loop over the sides to compute the one-sided slopes

  for (is = 0; is < nside; is++)
  {
    in = is + 1;

    // when the current element is an internal cell
    if (elem->neighbor(is) != NULL)
    {
      const Elem * neig = elem->neighbor(is);

      // _dof_map.dof_indices(neig, neig_dof_indices);
      const std::vector<dof_id_type> & neig_dof_indices = _dof_indices[neig->id()];

      // if (this->hasBlocks(neig->subdomain_id()))
      {
        xc[in] = _elem_centroid[neig->id()](0);

        // get the cell-average variable in this neighbor cell
        ucell[in][0] = (*_sys.solution)(neig_dof_indices[0]);

        // calculate the one-sided slopes of primitive variables

        for (iv = 0; iv < nvars; iv++)
          sigma[in][iv] = (ucell[0][iv] - ucell[in][iv]) / (xc[0] - xc[in]);
      }
      // else
      //   bflag = in;
    }
    // when the current element is at the boundary,
    // we choose not to construct the slope in 1D just for convenience.
    else
    {
      bflag = in;
    }
  }

  return ugrad;
}

void
RDGSystem::slopeReconstruction()
{
  // std::cerr << "RDGSystem::slopeReconstruction()" << std::endl;

  const MeshBase & mesh = this->mesh().getMesh();


  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;

    // TODO: slope reconstruction

    // slope limiting
    dof_id_type elem_id = elem->id();
    _lslope[elem_id] = limitElementSlope(elem);
    // std::cerr << "e=" << elem_id << ", slope = " << _lslope[elem_id][0] << std::endl;
  }
}

void
RDGSystem::assembleRHS()
{
  // std::cerr << "RDGSystem::assembleRHS()" << std::endl;

  const MeshBase & mesh = this->mesh().getMesh();
  // const unsigned int dim = mesh.mesh_dimension();

  FEType fe_type(CONSTANT, MONOMIAL);

  // UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  // UniquePtr<QBase> qrule(QBase::build(QGAUSS, dim, FIRST));
  // fe->attach_quadrature_rule(qrule.get());

  // UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  // UniquePtr<QBase>  qface (QBase::build(QGAUSS, dim - 1, FIRST));
  // fe_face->attach_quadrature_rule (qface.get());

  // const std::vector<Real> & JxW = fe->get_JxW();
  // const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // const std::vector<Real> & JxW = fe_face->get_JxW();
  // const std::vector<Point> & qface_point = fe_face->get_xyz();
  // const std::vector<Point> & normals = fe_face->get_normals();

  DenseVector<Number> rhs;
  rhs.resize(2);

  // std::vector<dof_id_type> dof_indices;
  // std::vector<dof_id_type> neig_dof_indices;

  std::vector<dof_id_type> dofs;
  dofs.resize(2);

  auto el_end = mesh.active_local_elements_end();
  for (auto el = mesh.active_local_elements_begin(); el != el_end; ++el)
  {
    const Elem * elem = *el;
    const dof_id_type elem_id = elem->id();

    // fe->reinit(elem);

    Point centroid = _elem_centroid[elem_id];

    // _dof_map.dof_indices(elem, dof_indices);
    const std::vector<dof_id_type> & dof_indices = _dof_indices[elem_id];
    dofs[0] = dof_indices[0];

    Real u_old = _sys.old_solution(dof_indices[0]);

    rhs(0) = _elem_JxW[elem_id] * (u_old / _dt);
    // std::cerr << "u_old / dt = " << rhs(0) << std::endl;
    _sys.rhs->add_vector(rhs, dof_indices);


    unsigned int nvars = 1;

    std::vector<RealGradient> ugrad(nvars, RealGradient(0., 0., 0.));
    RealGradient dvec;
    Real uc = (*_sys.solution)(dof_indices[0]);
    ugrad = _lslope[elem_id];

    for (unsigned int side = 0; side<elem->n_sides(); side++)
    {
      std::vector<BoundaryID> boundary_ids = _mesh.getBoundaryIDs(elem, side);

      // fe_face->reinit(elem, side);

      std::pair<dof_id_type, unsigned int> key(elem_id, side);

      const Point & dwave = _normals[key];
      const Point & face_point = _face_point[key];

      // unsigned int nvars = 1;

      // std::vector<RealGradient> ugrad(nvars, RealGradient(0., 0., 0.));
      RealGradient dvec;

      dvec = face_point - centroid;
      // calculate the variable at face center
      Real uvec1 = uc + ugrad[0] * dvec;

      if (boundary_ids.size() > 0)
        for (std::vector<BoundaryID>::iterator it = boundary_ids.begin(); it != boundary_ids.end(); ++it)
        {
          // onBoundary(elem, side, *it);

          // unsigned int nvars = 1;
          //
          // std::vector<RealGradient> ugrad(nvars, RealGradient(0., 0., 0.));
          // RealGradient dvec;
          // Real uc, un;

          // element
          // ugrad = _lslope[elem_id];
          // dvec = qface_point[0] - elem->centroid();
          // // calculate the variable at face center
          // uc = (*_sys.solution)(dof_indices[0]);
          // Real uvec1 = uc + ugrad[0] * dvec;

          // {
          //   const Point & dwave = normals[0];

            // assume the velocity vector is constant, e.g. = (1., 1., 1.)
            RealVectorValue uadv1(1.0, 1.0, 1.0);

            // finally calculate the flux
            Real flux = (uadv1 * dwave) * uvec1;

            rhs(0) = -flux * 1;

            // std::cerr << " BF: e=" << elem_id << ", (" << side <<"), u1=(" << uvec1 <<"), flux=(" << flux << "), normal=(" << dwave(0) << ")\n";

            _sys.rhs->add_vector(rhs, dof_indices);

          // }
        }

      if (elem->neighbor(side) != NULL)
      {
        // Pointer to the neighbor we are currently working on.
        const Elem * neighbor = elem->neighbor(side);

        // Get the global id of the element and the neighbor
        const dof_id_type neighbor_id = neighbor->id();

        // if ((neighbor->active() && (neighbor->level() == elem->level()) && (elem_id < neighbor_id)) || (neighbor->level() < elem->level()))
        if (elem_id < neighbor_id)
        {
          // _dof_map.dof_indices(neighbor, neig_dof_indices);
          const std::vector<dof_id_type> & neig_dof_indices = _dof_indices[neighbor_id];
          dofs[1] = neig_dof_indices[0];

          // unsigned int nvars = 1;

          // std::vector<RealGradient> ugrad(nvars, RealGradient(0., 0., 0.));
          // RealGradient dvec;
          Real un;

          // // element
          // ugrad = _lslope[elem_id];
          // dvec = qface_point[0] - elem->centroid();
          // // calculate the variable at face center
          // uc = (*_sys.solution)(dof_indices[0]);
          // Real uvec1 = uc + ugrad[0] * dvec;

          // neighbor
          ugrad = _lslope[neighbor_id];
          dvec = face_point - _elem_centroid[neighbor->id()];
          un = (*_sys.solution)(neig_dof_indices[0]);
          Real uvec2 = un + ugrad[0] * dvec;


          // const Point & dwave = normals[0];

          // assume a constant velocity on the left
          RealVectorValue uadv1(1.0, 1.0, 1.0);
          // assume a constant velocity on the right
          RealVectorValue uadv2(1.0, 1.0, 1.0);

          // normal velocity on the left and right
          Real vdon1 = uadv1 * dwave;
          Real vdon2 = uadv2 * dwave;

          // calculate the so-called a^plus and a^minus
          Real aplus = 0.5 * (vdon1 + std::abs(vdon1));
          Real amins = 0.5 * (vdon2 - std::abs(vdon2));

          // finally calculate the flux
          Real flux = aplus * uvec1 + amins * uvec2;

          rhs(0) = -flux * 1;
          rhs(1) =  flux * 1;

          // std::cerr << "ISF: e=" << elem_id << ", (" << side <<"), n=" << neighbor_id << ", u1=(" << uvec1 <<"), u2=(" << uvec2 << "), flux=(" << flux << ")\n";

          _sys.rhs->add_vector(rhs, dofs);
        }
      }
    } // sides
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

  _sys.rhs->zero();
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
  _sys.rhs->close();

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

  // sys().rhs->zero();
  // sys().rhs->close();
}
