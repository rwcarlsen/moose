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

/**
 * Explicit system to be solved
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

  virtual TransientExplicitSystem & sys() { return _sys; }

protected:
  class RDGAssembly : public System::Assembly
  {
  public:
    RDGAssembly(RDGSystem & rdg_system);

    virtual void assemble () override;

  protected:
    RDGSystem & _rdg_system;
  };

  RDGAssembly _rdg_assembly;
  TransientExplicitSystem & _sys;

  NumericVector<Number> & _mass_matrix;
};

#endif /* RDGSYSTEM_H */
