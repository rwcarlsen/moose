/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGSYSTEM_H
#define RDGSYSTEM_H

#include <unordered_map>

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

  struct elemprops {
    Real JxW;
    Point centroid;
    std::vector<dof_id_type> dof_indices;
    std::vector<Point> face_point;
    std::vector<Point> normals;
    std::vector<RealGradient> lslope;
  };
  std::vector<elemprops> _elemprops;
  elemprops& insertprop(dof_id_type elem_id)
  {
    if (_elemprops.size() <= elem_id)
      _elemprops.resize(elem_id+1, elemprops{});
    return _elemprops[elem_id];
  }

  DenseVector<Number> _rhs;
  std::vector<dof_id_type> _dofs;
  Point _centroid;
  RealGradient _dvec;

  RealVectorValue _uadv1;
  RealVectorValue _uadv2;


  friend class RDGProblem;
};

#endif /* RDGSYSTEM_H */
