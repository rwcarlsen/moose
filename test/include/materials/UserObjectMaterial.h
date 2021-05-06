//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class MaterialPropertyUserObject;

class UserObjectMaterial : public Material
{
public:
  static InputParameters validParams();

  UserObjectMaterial(const InputParameters & parameters);
  virtual void computeQpProperties();

protected:
  MaterialProperty<Real> & _prop;
  const MaterialPropertyUserObject & _user_object;
};
