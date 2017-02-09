/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef BOUNDARYFLUXBASE_H
#define BOUNDARYFLUXBASE_H

#include "SideUserObject.h"

// Forward Declarations
class BoundaryFluxBase;

template<>
InputParameters validParams<BoundaryFluxBase>();

/**
 * A base class for computing/caching fluxes at boundaries
 *
 * Notes:
 *
 *   1. When systems of equations are being solved, the fluxes are treated as vectors.
 *      To avoid recomputing the flux at the boundary, we compute it just once
 *      and then when it is needed, we just return the cached value.
 *
 *   2. Derived classes need to override `calcFlux` and `calcJacobian`.
 */
class BoundaryFluxBase : public SideUserObject
{
public:
  BoundaryFluxBase(const InputParameters & parameters);

  virtual void execute();
  virtual void initialize();
  virtual void finalize();
  virtual void threadJoin(const UserObject & uo);

  virtual void computeFlux();
  virtual void computeJacobian();

  virtual const std::vector<Real> & getFlux() const;
  /**
   * Get the boundary flux vector
   * @param[in]   iside     local  index of current side
   * @param[in]   ielem     global index of the current element
   * @param[in]   uvec1     vector of variables on the host side
   * @param[in]   dwave     vector of unit normal
   */
  // virtual const std::vector<Real> & getFlux(unsigned int iside,
  //                                           dof_id_type ielem,
  //                                           const std::vector<Real> & uvec1,
  //                                           const RealVectorValue & dwave,
  //                                           THREAD_ID tid) const;

  /**
   * Solve the Riemann problem on the boundary face
   * @param[in]   iside     local  index of current side
   * @param[in]   ielem     global index of the current element
   * @param[in]   uvec1     vector of variables on the host side
   * @param[in]   dwave     vector of unit normal
   * @param[out]  flux      flux vector for conservation equations
   */
  virtual void calcFlux(unsigned int iside,
                        dof_id_type ielem,
                        const std::vector<Real> & uvec1,
                        const RealVectorValue & dwave,
                        std::vector<Real> & flux) = 0;

  virtual const DenseMatrix<Real> & getJacobian() const;
  /**
   * Get the boundary Jacobian matrix
   * @param[in]   iside     local  index of current side
   * @param[in]   ielem     global index of the current element
   * @param[in]   uvec1     vector of variables on the host side
   * @param[in]   dwave     vector of unit normal
   */
  // virtual const DenseMatrix<Real> & getJacobian(unsigned int iside,
  //                                               dof_id_type ielem,
  //                                               const std::vector<Real> & uvec1,
  //                                               const RealVectorValue & dwave,
  //                                               THREAD_ID tid) const;

  /**
   * Compute the Jacobian matrix on the boundary face
   * @param[in]   iside     local  index of current side
   * @param[in]   ielem     global index of the current element
   * @param[in]   uvec1     vector of variables on the host side
   * @param[in]   dwave     vector of unit normal
   * @param[out]  jac1      Jacobian matrix contribution
   */
  virtual void calcJacobian(unsigned int iside,
                            dof_id_type ielem,
                            const std::vector<Real> & uvec1,
                            const RealVectorValue & dwave,
                            DenseMatrix<Real> & jac1) = 0;

protected:
  /// Threaded storage for fluxes
  std::vector<Real> _flux;

  /// Threaded storage for jacobians
  DenseMatrix<Real> _jac1;

private:
  static Threads::spin_mutex _mutex;
};

#endif // BOUNDARYFLUXBASE_H
