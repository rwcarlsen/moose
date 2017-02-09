/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGFLUXBC_H
#define RDGFLUXBC_H

#include "IntegratedBC.h"
#include "BoundaryFluxBase.h"

//Forward Declarations
class RDGFluxBC;

template<>
InputParameters validParams<RDGFluxBC>();

/**
 * A boundary condition that picks a component from a boundary flux
 */
class RDGFluxBC : public IntegratedBC
{
public:
  RDGFluxBC(const InputParameters & parameters);
  virtual ~RDGFluxBC() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  /// the component of the flux to use
  unsigned int _component;

  /// boundary flux object
  const BoundaryFluxBase & _flux;
};

#endif
