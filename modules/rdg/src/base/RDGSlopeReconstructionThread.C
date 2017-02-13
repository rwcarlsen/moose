/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RDGSlopeReconstructionThread.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "SwapBackSentinel.h"
#include "SlopeReconstructionBase.h"
#include "SlopeLimitingBase.h"

// libmesh includes
#include "libmesh/threads.h"

RDGSlopeReconstructionThread::RDGSlopeReconstructionThread(FEProblemBase & fe_problem,
                                                           const MooseObjectWarehouse<SlopeReconstructionBase> & sr_objects,
                                                           const MooseObjectWarehouse<SlopeLimitingBase> & limiting_objects)
  : ThreadedElementLoop<ConstElemRange>(fe_problem),
    _nl(fe_problem.getNonlinearSystemBase()),
    _sr_objects(sr_objects),
    _limiting_objects(limiting_objects)
{
}

// Splitting Constructor
RDGSlopeReconstructionThread::RDGSlopeReconstructionThread(RDGSlopeReconstructionThread & x, Threads::split split)
  : ThreadedElementLoop<ConstElemRange>(x, split),
    _nl(x._nl),
    _sr_objects(x._sr_objects),
    _limiting_objects(x._limiting_objects)
{
}

RDGSlopeReconstructionThread::~RDGSlopeReconstructionThread()
{
}

void
RDGSlopeReconstructionThread::subdomainChanged()
{
  _fe_problem.subdomainSetup(_subdomain, _tid);

  std::set<MooseVariable *> needed_moose_vars;
  _sr_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);

  _sr_objects.subdomainSetup(_subdomain, _tid);
  _limiting_objects.subdomainSetup(_subdomain, _tid);

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _limiting_objects.updateBlockVariableDependency(_subdomain, needed_moose_vars, _tid);

  // _fe_problem.prepareMaterials(_subdomain, _tid);
}

void
RDGSlopeReconstructionThread::onElement(const Elem * elem)
{
  _fe_problem.prepare(elem, _tid);
  // _fe_problem.reinitElem(elem, _tid);

  if (_sr_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<MooseSharedPointer<SlopeReconstructionBase> > & objects = _sr_objects.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & uo : objects)
      uo->reconstructElementSlope();
  }

  // slope limiting
  if (_limiting_objects.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<MooseSharedPointer<SlopeLimitingBase> > & objects = _limiting_objects.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & uo : objects)
      uo->computeSlopeLimiter();
  }
}

void
RDGSlopeReconstructionThread::onBoundary(const Elem *elem, unsigned int side, BoundaryID bnd_id)
{
}

void
RDGSlopeReconstructionThread::onInterface(const Elem *elem, unsigned int side, BoundaryID bnd_id)
{
}

void
RDGSlopeReconstructionThread::onInternalSide(const Elem *elem, unsigned int side)
{
}

void
RDGSlopeReconstructionThread::postElement(const Elem * /*elem*/)
{
}

void
RDGSlopeReconstructionThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
}

void
RDGSlopeReconstructionThread::join(const RDGSlopeReconstructionThread & /*y*/)
{
}
