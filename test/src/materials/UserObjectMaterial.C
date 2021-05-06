//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UserObjectMaterial.h"
#include "MaterialPropertyUserObject.h"

registerMooseObject("MooseTestApp", UserObjectMaterial);

InputParameters
UserObjectMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<UserObjectName>("user_object", "element integral user object to get value from");
  params.addParam<MaterialPropertyName>("mat_prop", "material property name");
  return params;
}

UserObjectMaterial::UserObjectMaterial(const InputParameters & parameters)
  : Material(parameters),
    _prop(declareProperty<Real>(getParam<MaterialPropertyName>("mat_prop"))),
    _user_object(getUserObject<MaterialPropertyUserObject>("user_object"))
{
  setRandomResetFrequency(EXEC_TIMESTEP_END);
}

void
UserObjectMaterial::computeQpProperties()
{
  _prop[_qp] = _user_object.getElementalValue(_current_elem->id());
}
