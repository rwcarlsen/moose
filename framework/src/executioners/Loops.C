
#include "Loops.h"

#include "NonlinearSystem.h"

template<>
InputParameters validParams<Loops>()
{
  InputParameters params = validParams<MooseObject>();
  params.registerBase("Loops");
  params.addParam<std::string>("flavor", "steady", "problem type (e.g. steady, transient, etc.)");
  return params;
}

Loops::Loops(const InputParameters & parameters) :
    MooseObject(parameters),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(parameters, "Executioners"),
    _problem(nullptr),
    _flavor(parameters.get<std::string>("flavor")),
    _max_time_steps(1)
{
  if (_flavor == "steady") _max_time_steps = 1;
}

void
Loops::initialize(QueenOfHearts& queen, FEProblem* problem)
{
  _problem = problem;

  LOOP_FROM_METHODS(queen, setup, teardown);
#ifdef LIBMESH_ENABLE_AMR
  LOOP_FROM_METHODS(queen, meshRefinementBegin, meshRefinementEnd);
#endif
  LOOP_FROM_METHODS(queen, timeStepBegin, timeStepEnd);
  LOOP_FROM_METHODS(queen, solveBegin, solveEnd);
}

bool Loops::setup(IterInfo info)
{
  if (_flavor == "steady" && _app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return true;
  }
  if (_flavor == "steady" && _problem->getNonlinearSystem().containsTimeKernel())
    mooseError("time kernels are not allowed in steady state simulations");

  _problem->initialSetup();
  _problem->outputStep(EXEC_INITIAL);

  if (!_app.isRecovering())
    _problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _problem.timeStep() = 1;
  // need to keep _time in sync with _time_step to get correct output
  _problem.time() = 1;

  return false;
}

bool Loops::teardown(IterInfo info) { return true;}

bool Loops::meshRefinementBegin(IterInfo info)
{
  return false;
}

bool Loops::meshRefinementEnd(IterInfo info)
{
  unsigned int max_steps = _problem->adaptivity().getSteps();
  if (info.curr_loop_iter < max_steps) _problem->adaptMesh();
  return info.curr_loop_iter >= max_steps;
}

bool Loops::timeStepBegin(IterInfo info)
{
  _problem->timestepSetup();
  _problem->execute(EXEC_TIMESTEP_BEGIN);
  _problem->outputStep(EXEC_TIMESTEP_BEGIN);

  // Update warehouse active objects
  _problem->updateActiveObjects();
  return false;
}

bool Loops::timeStepEnd(IterInfo info)
{
  if (_flavor == "steady" && !_problem->converged()) {
    _console << "aborting early: solve did not converge\n";
    return true;
  }

  _problem->onTimestepEnd();
  _problem->execute(EXEC_TIMESTEP_END);

  _problem->computeIndicators();
  _problem->computeMarkers();

  _problem->outputStep(EXEC_TIMESTEP_END);
  return info.curr_loop_iter >= _max_time_steps;
}

bool Loops::solveBegin(IterInfo info)
{
    _problem->solve();
    return false;
}

bool Loops::solveEnd(IterInfo info)
{
  return true;
}

