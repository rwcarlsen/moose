/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGSYSTEM_H
#define RDGSYSTEM_H

#include "SystemBase.h"
#include "NonlinearSystemBase.h"
#include "SlopeReconstructionBase.h"
#include "SlopeLimitingBase.h"
#include "BoundaryFluxBase.h"
#include "InternalSideFluxBase.h"
#include "libmesh/linear_implicit_system.h"


/**
 * A system that hold a RDG system being solved
 *
 */
class RDGSystem : public NonlinearSystemBase
{
public:
  RDGSystem(FEProblemBase & problem, const std::string & name);
  virtual ~RDGSystem();

  virtual void solve() override;

  /**
   * Quit the current solve as soon as possible.
   */
  virtual void stopSolve() override;

  /**
   * Returns the current nonlinear iteration number.  In libmesh, this is
   * updated during the nonlinear solve, so it should be up-to-date.
   */
  virtual unsigned int getCurrentNonlinearIterationNumber() override { return 0; }

  virtual void setupFiniteDifferencedPreconditioner() override;

  /**
   * Returns the convergence state
   * @return true if converged, otherwise false
   */
  virtual bool converged() override;

  virtual NumericVector<Number> & RHS() override { return *_sys.rhs; }

  virtual NonlinearSolver<Number> * nonlinearSolver() override { return NULL; }

  virtual NumericVector<Number> & solutionOld() override { return *_sys.old_local_solution; }

  virtual NumericVector<Number> & solutionOlder() override { return *_sys.older_local_solution; }

  virtual TransientLinearImplicitSystem & sys() { return _sys; }

  virtual void assemble();

  void needMatrix(bool state) { _need_matrix = state; }

protected:
  class RDGAssembly : public System::Assembly
  {
  public:
    RDGAssembly(RDGSystem & rdg_system);

    virtual void assemble() override;

  protected:
    RDGSystem & _rdg_system;
  };

  RDGAssembly _rdg_assembly;
  TransientLinearImplicitSystem & _sys;
  const DofMap & _dof_map;
  Real & _dt;
  bool _need_matrix;

  void assembleMassMatrix();
  void slopeReconstruction();
  void assembleRHS();

  std::vector<RealGradient> limitElementSlope(const Elem * elem);

  MooseObjectWarehouse<SlopeReconstructionBase> _reconstruction_objects;
  MooseObjectWarehouse<SlopeLimitingBase> _limiting_objects;
  MooseObjectWarehouse<BoundaryFluxBase> _boundary_flux_objects;
  MooseObjectWarehouse<InternalSideFluxBase> _internal_side_flux_objects;

  /// Stored JxW per element
  std::map<dof_id_type, Real> _elem_JxW;
  /// Stored element centroids
  std::map<dof_id_type, Point> _elem_centroid;
  std::map<dof_id_type, std::vector<dof_id_type> > _dof_indices;

  std::map<std::pair<dof_id_type, unsigned int>, Point> _face_point;
  std::map<std::pair<dof_id_type, unsigned int>, Point> _normals;

  /// store the updated slopes into this map indexed by element ID
  std::map<dof_id_type, std::vector<RealGradient> > _lslope;

  friend class RDGProblem;
};

#endif /* RDGSYSTEM_H */
