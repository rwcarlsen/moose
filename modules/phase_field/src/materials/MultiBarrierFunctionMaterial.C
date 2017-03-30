/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "MultiBarrierFunctionMaterial.h"

template <>
InputParameters
validParams<MultiBarrierFunctionMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Double well phase transformation barrier free energy contribution.");
  params.addParam<std::string>("function_name", "g", "actual name for g(eta_i)");
  MooseEnum h_order("SIMPLE=0", "SIMPLE");
  params.addParam<MooseEnum>(
      "g_order", h_order, "Polynomial order of the switching function h(eta)");
  params.addParam<bool>("well_only",
                        false,
                        "Make the g zero in [0:1] so it only contributes to "
                        "enforcing the eta range and not to the phase "
                        "transformation berrier.");
  params.addRequiredCoupledVar("etas", "eta_i order parameters, one for each h");
  return params;
}

MultiBarrierFunctionMaterial::MultiBarrierFunctionMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _function_name(getParam<std::string>("function_name")),
    _g_order(getParam<MooseEnum>("g_order")),
    _well_only(getParam<bool>("well_only")),
    _num_eta(coupledComponents("etas")),
    _eta(_num_eta),
{
  // declare derivative properties, fetch eta values
  bind_mat_prop(_function_name, computeProp);

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    auto eta_name = getVar("etas", i)->name();
    addPropFunc(derivProp(_function_name, eta_name),
                [this, i](const Location & loc) { return computePropd1(loc, i); });
    addPropFunc(derivProp(_function_name, eta_name, eta_name),
                [this, i](const Location & loc) { return computePropd2(loc, i); });
  }
}

Real
MultiBarrierFunctionMaterial::computeProp(const Location & loc)
{
  Real g = 0.0;
  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    const Real n = (*_eta[i])[loc.qp()];
    if (!(_well_only && n >= 0.0 && n <= 1.0) && _g_order == 0)
      g += n * n * (1.0 - n) * (1.0 - n);
  }
  return g;
}

Real
MultiBarrierFunctionMaterial::computePropd1(const Location & loc, unsigned int i)
{
  const Real n = (*_eta[i])[loc.qp()];
  if (_g_order != 0 || (_well_only && n >= 0.0 && n <= 1.0))
    return 0;
  return 2.0 * n * (n - 1.0) * (2.0 * n - 1.0);
}

Real
MultiBarrierFunctionMaterial::computePropd2(const Location & loc, unsigned int i)
{
  const Real n = (*_eta[i])[loc.qp()];
  if (_g_order != 0 || (_well_only && n >= 0.0 && n <= 1.0))
    return 0;
  return 12.0 * (n * n - n) + 2.0;
}
