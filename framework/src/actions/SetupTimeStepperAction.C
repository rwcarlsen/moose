/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SetupTimeStepperAction.h"
#include "Transient.h"
#include "Factory.h"
#include "TimeStepper.h"

template<>
InputParameters validParams<SetupTimeStepperAction>()
{
  InputParameters params = validParams<MooseObjectAction>();

  return params;
}

SetupTimeStepperAction::SetupTimeStepperAction(InputParameters parameters) :
    MooseObjectAction(parameters),
    _already_ran(false)
{
}

void
SetupTimeStepperAction::act()
{
  if (_already_ran) return;
  _already_ran = true;

  if (_problem->isTransient())
  {
    Transient * transient = dynamic_cast<Transient *>(_executioner.get());
    if (transient == NULL)
      mooseError("You can setup time stepper only with executioners of transient type.");
    _moose_object_pars.set<FEProblem *>("_fe_problem") = _problem.get();
    _moose_object_pars.set<Real>("_dt_min") = transient->dtMin();
    _moose_object_pars.set<Real>("_dt_max") = transient->dtMax();
    _moose_object_pars.set<Real>("_end_time") = transient->endTime();
    _moose_object_pars.set<bool>("_verbose") = transient->verbose();
    _moose_object_pars.set<bool>("_timestep_tolerance") = transient->timestepTol();

    if (_app._use_loops)
      _app._time_stepper = _factory.create<TimeStepper>(_type, "TimeStepper", _moose_object_pars);
    else
    {
      MooseSharedPointer<TimeStepper> ts = _factory.create<TimeStepper>(_type, "TimeStepper", _moose_object_pars);
      transient->setTimeStepper(ts);
    }
  }
}
