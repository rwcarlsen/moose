
#include "ExampleConvection.h"
#include <cmath>

template<>
InputParameters validParams<ExampleConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("kappa", "bla");
  return params;
}

ExampleConvection::ExampleConvection(const InputParameters & parameters) :
  Kernel(parameters),
   _kappa(getParam<Real>("kappa"))
{}

Real ExampleConvection::computeQpResidual() {
  return _grad_test[_i][_qp] *
             (1 / sqrt(1 + pow(_grad_u[_qp].norm(), 2)) * _grad_u[_qp]) +
         _kappa * _test[_i][_qp] * _u[_qp];
}

Real ExampleConvection::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
  return 0;
}

