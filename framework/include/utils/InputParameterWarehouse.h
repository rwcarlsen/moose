//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INPUTPARAMETERWAREHOUSE_H
#define INPUTPARAMETERWAREHOUSE_H

#include <gtest/gtest.h>
#include "MooseObjectName.h"
#include "MooseObjectParameterName.h"
#include "MooseTypes.h"
#include "Factory.h"
#include "ControlOutput.h"

// Forward declarations
class InputParameters;

/**
 * Storage container for all InputParamter objects.
 *
 * This object is responsible for InputParameter objects, all MooseObjects should
 * contain a reference to the parameters object stored here.
 *
 * To avoid abuse, this warehouse is also designed to restrict the ability to change the parameter
 * to Control objects only.
 */
class InputParameterWarehouse
{
public:
  /**
   * Class constructor
   */
  InputParameterWarehouse();

  /**
   * Destruction
   */
  virtual ~InputParameterWarehouse() = default;

  ///@{
  /**
   * Return a const reference to the InputParameters for the named object
   * @param tag The tag of the object (e.g., 'Kernel')
   * @param name The name of the parameters object, including the tag (name only input) or
   * MooseObjectName object
   * @param tid The thread id
   * @return A const reference to the warehouse copy of the InputParameters
   */
  const InputParameters & getInputParametersObject(const std::string & name,
                                                   THREAD_ID tid = 0) const;
  const InputParameters & getInputParametersObject(const std::string & tag,
                                                   const std::string & name,
                                                   THREAD_ID tid = 0) const;
  const InputParameters & getInputParametersObject(const MooseObjectName & object_name,
                                                   THREAD_ID tid = 0) const;
  ///@}
  /**
   * Return const reference to the map containing the InputParameter objects
   */
  const std::multimap<MooseObjectName, std::shared_ptr<InputParameters>> &
  getInputParameters(THREAD_ID tid = 0) const;

  /***
   * Helper method for printing controllable items.
   */
  std::string dumpChangedControls(bool reset_changed) const;

  /**
   * Returns a copy of the current values for a controllable parameter. This method
   * is designed to provide access to objects for monitoring the state of a controllable parameter.
   */
  template <typename T>
  std::vector<T> getControllableParameterValues(const MooseObjectParameterName & input) const;

private:
  /// Storage for the InputParameters objects
  /// TODO: Remove multimap
  std::vector<std::multimap<MooseObjectName, std::shared_ptr<InputParameters>>> _input_parameters;

  /**
   * Method for adding a new InputParameters object
   * @param parameters The InputParameters object to copy and store in the warehouse
   * @return A reference to the warehouse copy of the InputParameters, this
   *         is what should be passed into the MooseObjects constructors.
   *
   * A new object is created from the old object because InputParameters objects
   * are generic until Factory::create() is called and the actual MooseObject
   * is created.
   *
   * This method is private, because only the factories that are creating objects should be
   * able to call this method.
   */
  InputParameters &
  addInputParameters(const std::string & name, InputParameters parameters, THREAD_ID tid = 0);

  ///@{
  /**
   * Return a reference to the InputParameters for the named object
   * @param tag The tag of the object (e.g., 'Kernel')
   * @param name The name of the parameters object, including the tag (name only input) or
   * MooseObjectName object
   * @param tid The thread id
   * @return A reference to the warehouse copy of the InputParameters
   *
   * If you are using this method to access a writable reference to input parameters, this
   * will break the ability to control the parameters with the MOOSE control logic system.
   * Only change parameters if you know what you are doing. Hence, this is private for a reason.
   */
  InputParameters & getInputParameters(const std::string & name, THREAD_ID tid = 0) const;
  InputParameters &
  getInputParameters(const std::string & tag, const std::string & name, THREAD_ID tid = 0) const;
  InputParameters & getInputParameters(const MooseObjectName & object_name,
                                       THREAD_ID tid = 0) const;
  ///@{

  /// The factory is allowed to call addInputParameters.
  friend MooseObjectPtr
  Factory::create(const std::string &, const std::string &, InputParameters, THREAD_ID, bool);

  /// Only controls are allowed to call getControllableParameter. The
  /// Control::getControllableParameter is the only method that calls getControllableParameter.
  /// However, this method cannot be made a friend explicitly because the method would need to be
  /// public. If the method is public then it is possible to call the method by getting access to
  /// the Control object.
  friend class Control;

};

#endif // INPUTPARAMETERWAREHOUSE_H
