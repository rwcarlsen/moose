
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
    _root(nullptr),
    _flavor(parameters.get<std::string>("flavor"))
{
}

void
Loops::initialize(FEProblem* problem)
{
  _problem = problem;
  ExecLoop* curr = nullptr;

  {
    SetupLoop* loop = new SetupLoop();
    loop->_steady = _flavor == "steady";
    curr = loop;
    _root = loop;
  }

#ifdef LIBMESH_ENABLE_AMR
  {
    MeshRefinementLoop* loop = new MeshRefinementLoop();
    curr->addChild(loop);
    curr = loop;
  }
#endif

  {
    TimeLoop* loop = new TimeLoop();
    loop->_steady = _flavor == "steady";
    curr->addChild(loop);
    curr = loop;
  }

  {
    SolveLoop* loop = new SolveLoop();
    curr->addChild(loop);
    curr = loop;
  }
}

void
Loops::run() {
  LoopContext ctx = {&_app, _problem};
  _root->run(ctx);
}

////////////////// SetupLoop ////////////////////
SetupLoop::SetupLoop() : _steady(true) { }

std::string SetupLoop::name()
{
  return "setup";
}

bool SetupLoop::beginIter(LoopContext& ctx)
{
  if (_steady && ctx.app->isRecovering())
  {
    //ctx.console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return true;
  }
  if (_steady && ctx.problem->getNonlinearSystem().containsTimeKernel())
    mooseError("time kernels are not allowed in steady state simulations");

  ctx.problem->initialSetup();
  ctx.problem->outputStep(EXEC_INITIAL);

  if (!ctx.app->isRecovering())
    ctx.problem->advanceState();

  if (_steady) {
    // first step in any steady state solve is always 1 (preserving backwards compatibility)
    ctx.problem->timeStep() = 1;
    // need to keep _time in sync with _time_step to get correct output
    ctx.problem->time() = 1;
  }

  return false;
}

bool SetupLoop::endIter(LoopContext& ctx)
{
  return true;
}

////////////////// MeshRefinementLoop ////////////////////
std::string MeshRefinementLoop::name()
{
  return "mesh-refinement-loop";
}

bool MeshRefinementLoop::beginIter(LoopContext& ctx)
{
  return false;
}

bool MeshRefinementLoop::endIter(LoopContext& ctx)
{
  unsigned int max_steps = ctx.problem->adaptivity().getSteps();
  if (iter() < max_steps) ctx.problem->adaptMesh();
  return iter() >= max_steps;
}

////////////////// TimeLoop ////////////////////
std::string TimeLoop::name()
{
  return "time-loop";
}

TimeLoop::TimeLoop() : _steady(true) { }

bool TimeLoop::beginIter(LoopContext& ctx)
{
  ctx.problem->timestepSetup();
  ctx.problem->execute(EXEC_TIMESTEP_BEGIN);
  ctx.problem->outputStep(EXEC_TIMESTEP_BEGIN);

  // Update warehouse active objects
  ctx.problem->updateActiveObjects();
  return false;
}

bool TimeLoop::endIter(LoopContext& ctx)
{
  if (_steady && !ctx.problem->converged()) {
    //ctx.console << "aborting early: solve did not converge\n";
    return true;
  }

  ctx.problem->onTimestepEnd();
  ctx.problem->execute(EXEC_TIMESTEP_END);

  ctx.problem->computeIndicators();
  ctx.problem->computeMarkers();

  ctx.problem->outputStep(EXEC_TIMESTEP_END);
  return _steady;
}

////////////////// SolveLoop ////////////////////
std::string SolveLoop::name()
{
  return "solve";
}

bool SolveLoop::beginIter(LoopContext& ctx)
{
    ctx.problem->solve();
    return false;
}

bool SolveLoop::endIter(LoopContext& ctx)
{
  return true;
}

