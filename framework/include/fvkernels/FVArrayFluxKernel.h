//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVKernel.h"
#include "FVUtils.h"
#include "NeighborCoupleable.h"
#include "TwoMaterialPropertyInterface.h"
#include "NeighborMooseVariableInterface.h"
#include "NeighborCoupleableMooseVariableDependencyIntermediateInterface.h"

class FaceInfo;

/// FVArrayFluxKernel is used for calculating residual contributions from numerical
/// fluxes from surface integral terms in a finite volume discretization of a
/// PDE (i.e.  terms where the divergence theorem is applied).  As with finite
/// element kernels, all solution values and material properties must be
/// indexed using the _qp member.  Note that all interfaces for finite volume
/// kernels are AD-based - be sure to use AD material properties and other AD
/// values to maintain good jacobian/derivative quality.
class FVArrayFluxKernel : public FVFluxKernelBase,
                          public NeighborMooseVariableInterface<RealEigenVector>
{
public:
  static InputParameters validParams();
  FVArrayFluxKernel(const InputParameters & params);

  /// Usually you should not override these functions - they have some super
  /// tricky stuff in them that you don't want to mess up!
  // @{
  virtual void computeResidual(const FaceInfo & fi) override;
  virtual void computeJacobian(const FaceInfo & fi) override;
  /// @}

protected:
  /// This is the primary function that must be implemented for flux kernel
  /// terms.  Material properties will be initialized on the face - using any
  /// reconstructed fv variable gradients if any.  Values for the solution are
  /// provided for both the elem and neighbor side of the face.
  virtual ADRealEigenVector computeQpResidual() = 0;

  /// Calculates and returns "grad_u dot normal" on the face to be used for
  /// diffusive terms.  If using any cross-diffusion corrections, etc. all
  /// those calculations will be handled for appropriately by this function.
  virtual ADRealEigenVector gradUDotNormal() const;

  virtual MooseVariableFieldBase & fieldVar() override { return _var; }

  MooseVariableFVArray & _var;

  /// The elem solution value of the kernel's _var for the current face.
  const ADArrayVariableValue & _u_elem;
  /// The neighbor solution value of the kernel's _var for the current face.
  const ADArrayVariableValue & _u_neighbor;
  /// The elem solution gradient of the kernel's _var for the current face.
  /// This is zero unless higher order reconstruction is used.
  const ADArrayVariableGradient & _grad_u_elem;
  /// The neighbor solution gradient of the kernel's _var for the current face.
  /// This is zero unless higher order reconstruction is used.
  const ADArrayVariableGradient & _grad_u_neighbor;
};
