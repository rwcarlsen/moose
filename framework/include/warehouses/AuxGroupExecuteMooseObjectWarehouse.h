//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef AUXGROUPEXECUTEMOOSEOBJECTWAREHOUSEBASE_H
#define AUXGROUPEXECUTEMOOSEOBJECTWAREHOUSEBASE_H

// MOOSE includes
#include "ExecuteMooseObjectWarehouse.h"

class UserObject;

/**
 * General warehouse for storing MooseObjects based on relation to IC and AuxKernel execution.
 */
template <typename T>
class AuxGroupExecuteMooseObjectWarehouse : public ExecuteMooseObjectWarehouse<T>
{

public:
  /// Using these from base class
  using MooseObjectWarehouse<T>::checkThreadID;
  using ExecuteMooseObjectWarehouse<T>::_all_objects;
  using ExecuteMooseObjectWarehouse<T>::_execute_objects;
  using ExecuteMooseObjectWarehouse<T>::_num_threads;

  /**
   * Constructor.
   */
  AuxGroupExecuteMooseObjectWarehouse(const ExecFlagEnum & flags, bool thread = true);

  /**
   * Access the AuxGroup via bracket operator.
   */
  const ExecuteMooseObjectWarehouse<T> & operator[](Moose::AuxGroup group) const;

  /**
   * Call this to separate the stored objects into the various AuxGroup categories.
   *
   * @see FEProblemBase::initialSetup()
   */
  void updateDependObjects(const std::set<std::string> & depend_ic,
                           const std::set<std::string> & depend_aux,
                           THREAD_ID tid = 0);

  /**
   * Performs a sort using the DependencyResolver.
   */
  void sort(THREAD_ID tid = 0);

  /**
   * Updates the various active lists of objects.
   */
  virtual void updateActive(THREAD_ID tid = 0) override;

protected:
  /// Storage for the group sorted objects (ALL is stored in the base class)
  std::vector<ExecuteMooseObjectWarehouse<T>> _group_objects;
};

template <typename T>
AuxGroupExecuteMooseObjectWarehouse<T>::AuxGroupExecuteMooseObjectWarehouse(
    const ExecFlagEnum & flags, bool threaded)
  : ExecuteMooseObjectWarehouse<T>(flags, threaded),
    _group_objects(3, ExecuteMooseObjectWarehouse<T>(flags, threaded)) // initialize group storage
{
}

template <typename T>
const ExecuteMooseObjectWarehouse<T> & AuxGroupExecuteMooseObjectWarehouse<T>::
operator[](Moose::AuxGroup group) const
{
  if (group == Moose::ALL)
    return *this;
  return _group_objects[group];
}

// determine when we need to run user objects based on whether any initial conditions or aux
// kernels depend on the user objects.  If so we need to run them either before ics, before aux
// kernels, or after aux kernels (if nothing depends on them).  Mark/store this information as
// attributes in the warehouse for later reference.
template <typename T>
void
groupUserObjects(TheWarehouse & w,
                 const std::vector<T *> & objs,
                 const std::set<std::string> & ic_deps,
                 const std::set<std::string> & aux_deps)
{
  // Notes about how this information is used later during the simulation:
  // We only need to run pre-ic objects for their "initial" exec flag time (not the others).
  //
  // For pre/post aux objects:
  //
  //     If an object was not run as a pre-ic object or it is a pre-ic object and
  //     shouldDuplicateInitialExecution returns true:
  //         * run the object at all its exec flag times.
  //     Else
  //         * run the object at all its exec flag times *except* "initial"
  //
  for (const auto obj : objs)
  {
    std::vector<Attribute> attribs;
    auto b = w.build();

    if (ic_deps.count(obj->name()) > 0)
      b.pre_ic(true);

    if ((obj->isParamValid("force_preaux") && obj->template getParam<bool>("force_preaux")) ||
        aux_deps.count(obj->name()) > 0 || ic_deps.count(obj->name()) > 0)
      b.pre_aux(true);
    else
      b.pre_aux(false);

    w.update(obj, b.attribs());
  }
}

template <typename T>
void
AuxGroupExecuteMooseObjectWarehouse<T>::updateDependObjects(
    const std::set<std::string> & depend_ic,
    const std::set<std::string> & depend_aux,
    THREAD_ID tid)
{
  checkThreadID(tid);

  const std::uint16_t initial_flag_mask = static_cast<std::uint16_t>(EXEC_INITIAL);
  const std::uint16_t not_initial_flag_mask = ~static_cast<std::uint16_t>(EXEC_INITIAL);
  const std::uint16_t all_flags = std::numeric_limits<std::uint16_t>::max();

  for (const auto & object_ptr : _all_objects[tid])
  {
    bool already_added = false;
    if (depend_ic.find(object_ptr->name()) != depend_ic.end())
    {
      _group_objects[Moose::PRE_IC].addObjectMask(object_ptr, tid, initial_flag_mask);
      already_added = !object_ptr->shouldDuplicateInitialExecution();
    }

    std::uint16_t remaining_flags = already_added ? not_initial_flag_mask : all_flags;
    if ((object_ptr->isParamValid("force_preaux") &&
         object_ptr->template getParam<bool>("force_preaux")) ||
        depend_aux.find(object_ptr->name()) != depend_aux.end() ||
        depend_ic.find(object_ptr->name()) != depend_ic.end())
      _group_objects[Moose::PRE_AUX].addObjectMask(object_ptr, tid, remaining_flags);
    else
      _group_objects[Moose::POST_AUX].addObjectMask(object_ptr, tid, remaining_flags);
  }
}

template <typename T>
void
AuxGroupExecuteMooseObjectWarehouse<T>::sort(THREAD_ID tid /*= 0*/)
{
  ExecuteMooseObjectWarehouse<T>::sort(tid);
  _group_objects[Moose::PRE_IC].sort(tid);
  _group_objects[Moose::PRE_AUX].sort(tid);
  _group_objects[Moose::POST_AUX].sort(tid);
}

template <typename T>
void
AuxGroupExecuteMooseObjectWarehouse<T>::updateActive(THREAD_ID tid /*=0*/)
{
  ExecuteMooseObjectWarehouse<T>::updateActive(tid);
  _group_objects[Moose::PRE_IC].updateActive(tid);
  _group_objects[Moose::PRE_AUX].updateActive(tid);
  _group_objects[Moose::POST_AUX].updateActive(tid);
}

#endif // AUXGROUPEXECUTEMOOSEOBJECTWAREHOUSEBASE_H
