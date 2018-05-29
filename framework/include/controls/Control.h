//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CONTROL_H
#define CONTROL_H

// MOOSE includes
#include "MooseObject.h"
#include "ControllableParameter.h"
#include "TransientInterface.h"
#include "SetupInterface.h"
#include "FunctionInterface.h"
#include "UserObjectInterface.h"
#include "PostprocessorInterface.h"
#include "VectorPostprocessorInterface.h"

// Forward declarations
class Control;
class FEProblemBase;
class InputParameterWarehouse;

template <>
InputParameters validParams<Control>();

/**
 * Base class for Control objects.
 *
 * These objects are create by the [Controls] block in the input file after
 * all other MooseObjects are created, so they have access to parameters
 * in all other MooseObjects.
 */
class Control : public MooseObject,
                public TransientInterface,
                public SetupInterface,
                public FunctionInterface,
                public UserObjectInterface,
                public Restartable,
                protected PostprocessorInterface,
                protected VectorPostprocessorInterface
{
public:
  /**
   * Class constructor
   * @param parameters The input parameters for this control object
   */
  Control(const InputParameters & parameters);

  /**
   * Class destructor
   */
  virtual ~Control() {}

  /**
   * Execute the control. This must be overridden.
   */
  virtual void execute() = 0;

  /**
   * (DEPRECATED) Return the valid "execute_on" options for Control objects
   */
  static MultiMooseEnum getExecuteOptions();

  /**
   * Return the Controls that must run before this Control
   */
  std::vector<std::string> & getDependencies() { return _depends_on; }

protected:
  template <typename T>
  void setControllableValue(const MooseObjectParameterName & name, const T & value);

private:
  /// A reference to the InputParameterWarehouse which is used for access the parameter objects
  InputParameterWarehouse & _input_parameter_warehouse;

};

template <typename T>
void
Control::setControllableValue(const MooseObjectParameterName & desired, const T & value)
{
  auto & p = _input_parameter_warehouse.getInputParameters(object_name);
  p.set<T>(object_name.parameter()) = value;
}

#endif
