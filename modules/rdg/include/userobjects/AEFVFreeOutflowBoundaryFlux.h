/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef AEFVFREEOUTFLOWBOUNDARYFLUX_H
#define AEFVFREEOUTFLOWBOUNDARYFLUX_H

#include "BoundaryFluxBase.h"

// Forward Declarations
class AEFVFreeOutflowBoundaryFlux;

template<>
InputParameters validParams<AEFVFreeOutflowBoundaryFlux>();

/**
 * Free outflow BC based boundary flux user object
 * for the advection equation
 * using a cell-centered finite volume method
 */
class AEFVFreeOutflowBoundaryFlux : public BoundaryFluxBase
{
public:

  AEFVFreeOutflowBoundaryFlux(const InputParameters & parameters);
  virtual ~AEFVFreeOutflowBoundaryFlux();

  virtual void computeFlux() override;
  virtual void computeJacobian() override;

  virtual void calcFlux(unsigned int iside,
                        dof_id_type ielem,
                        const std::vector<Real> & uvec1,
                        const RealVectorValue & dwave,
                        std::vector<Real> & flux) override;

  virtual void calcJacobian(unsigned int iside,
                            dof_id_type ielem,
                            const std::vector<Real> & uvec1,
                            const RealVectorValue & dwave,
                            DenseMatrix<Real> & jac1) override;

protected:
  // "1" denotes variable value from the host element

  /// piecewise constant variable values in host element
  const VariableValue &  _uc1;

  /// extrapolated variable values at side center
  const MaterialProperty<Real> &  _u1;
};

#endif
