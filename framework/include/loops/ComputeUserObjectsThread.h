//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEUSEROBJECTSTHREAD_H
#define COMPUTEUSEROBJECTSTHREAD_H

// MOOSE includes
#include "ThreadedElementLoop.h"

#include "libmesh/elem_range.h"

// libMesh forward declarations
namespace libMesh
{
template <typename T>
class NumericVector;
}

/**
 * Class for threaded computation of UserObjects.
 */
class ComputeUserObjectsThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeUserObjectsThread(FEProblemBase & problem, SystemBase & sys);
  // Splitting Constructor
  ComputeUserObjectsThread(ComputeUserObjectsThread & x, Threads::split);

  virtual ~ComputeUserObjectsThread();

  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void post() override;
  virtual void subdomainChanged() override;

  void join(const ComputeUserObjectsThread & /*y*/);

protected:
  const NumericVector<Number> & _soln;

private:
  template <typename T>
  TheWarehouse::Builder querySubdomain(Interfaces iface, std::vector<T> & results)
  {
    _fe_problem.theWarehouse()
        .build()
        .thread(_tid)
        .subdomain(_subdomain)
        .interfaces(iface)
        .queryInto(results);
  }
  template <typename T>
  TheWarehouse::Builder queryBoundary(Interfaces iface, BoundaryID bnd, std::vector<T> & results)
  {
    _fe_problem.theWarehouse().build().thread(_tid).boundary(bnd).interfaces(iface).queryInto(
        results);
  }
};

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

#endif // COMPUTEUSEROBJECTSTHREAD_H
