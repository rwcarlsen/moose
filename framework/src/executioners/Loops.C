
#include "Loops.h"

#include "NonlinearSystem.h"

void validParamsTimeLoop(InputParameters& params);

template<>
InputParameters validParams<Loops>()
{
  InputParameters params = validParams<MooseObject>();
  params.registerBase("Loops");
  params.addParam<std::string>("flavor", "steady-state", "problem type (e.g. steady, transient, etc.)");
  validParamsTimeLoop(params);
  return params;
}

Loops::Loops(const InputParameters & parameters) :
    MooseObject(parameters),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(parameters, "Executioners"),
    _problem(nullptr),
    _root(nullptr)
{
}

template <class T>
T* nestLoopUnder(ExecLoop* parent, const InputParameters& params, LoopContext* ctx) {
  T* loop = new T(params, ctx);
  if (parent != NULL)
    parent->addChild(loop);
  return loop;
}

void
Loops::initialize(FEProblem* problem)
{
  _problem = problem;
  LoopContext ctx = {&_app, _problem};
  ExecLoop* curr = nullptr;

  std::string flavor = _pars.get<std::string>("flavor");
  if (flavor == "steady-state")
  {
    curr = nestLoopUnder<SetupLoop>(NULL, _pars, &ctx);
    _root = curr;

#ifdef LIBMESH_ENABLE_AMR
    curr = nestLoopUnder<MeshRefinementLoop>(curr, _pars, &ctx);
#endif
    curr = nestLoopUnder<TimeLoop>(curr, _pars, &ctx);
    curr = nestLoopUnder<SolveLoop>(curr, _pars, &ctx);
  }
}

void
Loops::run() {
  LoopContext ctx = {&_app, _problem};
  _root->run(&ctx);
}

////////////////// SetupLoop ////////////////////
SetupLoop::SetupLoop(const InputParameters& params, LoopContext* ctx) : _steady(true)
{
  _steady = params.get<std::string>("flavor") == "steady-state";
}

std::string SetupLoop::name()
{
  return "setup";
}

bool SetupLoop::beginIter(LoopContext* ctx)
{
  if (_steady && ctx->app->isRecovering())
  {
    //ctx->console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return true;
  }
  if (_steady && ctx->problem->getNonlinearSystem().containsTimeKernel())
    mooseError("time kernels are not allowed in steady state simulations");

  ctx->problem->initialSetup();
  ctx->problem->outputStep(EXEC_INITIAL);

  if (!ctx->app->isRecovering())
    ctx->problem->advanceState();

  if (_steady)
  {
    // first step in any steady state solve is always 1 (preserving backwards compatibility)
    ctx->problem->timeStep() = 1;
    // need to keep _time in sync with _time_step to get correct output
    ctx->problem->time() = 1;
  }

  return false;
}

bool SetupLoop::endIter(LoopContext* ctx)
{
  return true;
}

////////////////// MeshRefinementLoop ////////////////////
MeshRefinementLoop::MeshRefinementLoop(const InputParameters& params, LoopContext* ctx) { }

std::string MeshRefinementLoop::name()
{
  return "mesh-refinement-loop";
}

bool MeshRefinementLoop::beginIter(LoopContext* ctx)
{
  return false;
}

bool MeshRefinementLoop::endIter(LoopContext* ctx)
{
  unsigned int max_steps = ctx->problem->adaptivity().getSteps();
  if (iter() < max_steps) ctx->problem->adaptMesh();
  return iter() >= max_steps;
}

void
validParamsTimeLoop(InputParameters& params)
{
  params.addParam<unsigned int>("num_steps",       std::numeric_limits<unsigned int>::max(),     "The number of timesteps in a transient run");
}

////////////////// TimeLoop ////////////////////
TimeLoop::TimeLoop(const InputParameters& params, LoopContext* ctx) :
    _num_steps(params.get<unsigned int>("num_steps")),
    _steps_taken(0),

    _t_step(ctx->problem->timeStep()),
    _time(ctx->problem->time()),
    _time_old(ctx->problem->timeOld()),
    _dt(ctx->problem->dt()),
    _dt_old(ctx->problem->dtOld()),
    _sync_times(ctx->app->getOutputWarehouse().getSyncTimes())
{
  if (params.get<std::string>("flavor") == "steady-state")
    _num_steps = 1;

  //////// from transient solver: ////////////
  //_t_step = 0;
  //_dt = 0;
  //_next_interval_output_time = 0.0;

  //// Either a start_time has been forced on us, or we want to tell the App about what our start time is (in case anyone else is interested.
  //if (_app.hasStartTime())
  //  _start_time = _app.getStartTime();
  //else if (parameters.isParamSetByUser("start_time"))
  //  _app.setStartTime(_start_time);

  //_time = _time_old = _start_time;
  //_problem.transient(true);

  //setupTimeIntegrator();

  //if (_app.halfTransient()) // Cut timesteps and end_time in half...
  //{
  //  _end_time /= 2.0;
  //  _num_steps /= 2.0;

  //  if (_num_steps == 0) // Always do one step in the first half
  //    _num_steps = 1;
  //}
}

std::string TimeLoop::name()
{
  return "time-loop";
}

bool TimeLoop::beginIter(LoopContext* ctx)
{
  ctx->problem->timestepSetup();
  ctx->problem->execute(EXEC_TIMESTEP_BEGIN);
  ctx->problem->outputStep(EXEC_TIMESTEP_BEGIN);

  // Update warehouse active objects
  ctx->problem->updateActiveObjects();
  return false;
}

bool TimeLoop::endIter(LoopContext* ctx)
{
  if (_num_steps == 1 && !ctx->problem->converged()) {
    //ctx->console << "aborting early: solve did not converge\n";
    return true;
  }

  ctx->problem->onTimestepEnd();
  ctx->problem->execute(EXEC_TIMESTEP_END);

  ctx->problem->computeIndicators();
  ctx->problem->computeMarkers();

  ctx->problem->outputStep(EXEC_TIMESTEP_END);
  return _num_steps >= iter();
}

////////////////// SolveLoop ////////////////////
std::string SolveLoop::name()
{
  return "solve";
}

bool SolveLoop::beginIter(LoopContext* ctx)
{
    ctx->problem->solve();
    return false;
}

bool SolveLoop::endIter(LoopContext* ctx)
{
  return true;
}

