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

#include "CreateLoopsAction.h"
#include "Factory.h"
#include "MooseApp.h"
#include "Loops.h"

template<>
InputParameters validParams<CreateLoopsAction>()
{
  InputParameters params = validParams<MooseObjectAction>();
  return params;
}


CreateLoopsAction::CreateLoopsAction(InputParameters params) :
    MooseObjectAction(params)
{
}

void
CreateLoopsAction::act()
{
  _moose_object_pars.set<FEProblem *>("_fe_problem") = _problem.get();
  _app._use_queen = true;
  _app._loops = _factory.create<Loops>(_type, "Loops", _moose_object_pars);
  _app._loops->initialize(_problem.get());
}
