/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGFLUX_H
#define RDGFLUX_H

#include "DGKernel.h"

class RDGFlux;
class InternalSideFluxBase;

template<>
InputParameters validParams<RDGFlux>();

/**
 * A DG kernel that picks a component of a flux
 */
class RDGFlux : public DGKernel
{
public:
  RDGFlux(const InputParameters & parameters);
  virtual ~RDGFlux();

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  /// The component of the flux
  unsigned int _component;

  /// flux user object
  const InternalSideFluxBase & _flux;
};

#endif
