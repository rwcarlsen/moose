/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef AEFVUpwindINTERNALSIDEFLUX_H
#define AEFVUpwindINTERNALSIDEFLUX_H

#include "InternalSideFluxBase.h"

// Forward Declarations
class AEFVUpwindInternalSideFlux;

template<>
InputParameters validParams<AEFVUpwindInternalSideFlux>();

/**
 * Upwind numerical flux scheme
 * for the advection equation
 * using a cell-centered finite volume method
 */
class AEFVUpwindInternalSideFlux : public InternalSideFluxBase
{
public:
  AEFVUpwindInternalSideFlux(const InputParameters & parameters);
  virtual ~AEFVUpwindInternalSideFlux();

  virtual void computeFlux() override;
  virtual void computeJacobian() override;

  virtual void calcFlux(unsigned int iside,
                        unsigned int ielem,
                        unsigned int ineig,
                        const std::vector<Real> & uvec1,
                        const std::vector<Real> & uvec2,
                        const RealVectorValue & dwave,
                        std::vector<Real> & flux) override;

  virtual void calcJacobian(unsigned int iside,
                            unsigned int ielem,
                            unsigned int ineig,
                            const std::vector<Real> & uvec1,
                            const std::vector<Real> & uvec2,
                            const RealVectorValue & dwave,
                            DenseMatrix<Real> & jac1,
                            DenseMatrix<Real> & jac2) override;

protected:
  // "1" denotes the "left" state
  // "2" denotes the "right" state

  /// piecewise constant variable values in cells
  const VariableValue &  _uc1;
  const VariableValue &  _uc2;

  /// extrapolated variable values at side center
  const MaterialProperty<Real> &  _u1;
  const MaterialProperty<Real> &  _u2;

};

#endif
