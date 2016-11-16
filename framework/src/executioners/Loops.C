
#include "Loops.h"

#include "NonlinearSystem.h"

void validParamsSetupLoop(InputParameters& params);
void validParamsSolveLoop(InputParameters& params);
void validParamsTimeLoop(InputParameters& params);
void validParamsPicardLoop(InputParameters& params);
void validParamsXfemLoop(InputParameters& params);

template<>
InputParameters validParams<Loops>()
{
  InputParameters params = validParams<MooseObject>();
  params.registerBase("Loops");
  params.addParam<std::string>("flavor", "steady-state", "problem type (e.g. steady, transient, etc.)");
  validParamsXfemLoop(params);
  validParamsSetupLoop(params);
  validParamsPicardLoop(params);
  validParamsSolveLoop(params);
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

void
Loops::initialize(FEProblem* problem)
{
  _problem = problem;
  LoopContext ctx(_app, *_problem);
  ExecLoop* curr = nullptr;

  std::string flavor = _pars.get<std::string>("flavor");
  if (flavor == "steady-state")
  {
    _root = new SetupLoop(_pars, ctx);
    curr = _root;

#ifdef LIBMESH_ENABLE_AMR
    curr = curr->addChild(new MeshRefinementLoop(_pars, ctx));
#endif
    curr = curr->addChild(new TimeLoop(_pars, ctx));
    curr = curr->addChild(new SolveLoop(_pars, ctx));
  }
}

void
Loops::run() {
  LoopContext ctx(_app, *_problem);
  _root->run(ctx);
}

////////////////// SetupLoop ////////////////////
void
validParamsSetupLoop(InputParameters& params)
{
}

SetupLoop::SetupLoop(const InputParameters& params, LoopContext& ctx) : _steady(true)
{
  _steady = params.get<std::string>("flavor") == "steady-state";
}

std::string SetupLoop::name()
{
  return "setup-loop";
}

void
SetupLoop::beginIter(LoopContext& ctx)
{
  std::cout << "setup begin\n";
  if (_steady && ctx.app().isRecovering())
    mooseError("cannot recover steady solves");
  if (_steady && ctx.problem().getNonlinearSystem().containsTimeKernel())
    mooseError("time kernels are not allowed in steady state simulations");

  ctx.problem().initialSetup();
  ctx.problem().outputStep(EXEC_INITIAL);

  if (!ctx.app().isRecovering())
    ctx.problem().advanceState();

  if (_steady)
  {
    // first step in any steady state solve is always 1 (preserving backwards compatibility)
    ctx.problem().timeStep() = 1;
    // need to keep _time in sync with _time_step to get correct output
    ctx.problem().time() = 1;
  }
}

void
SetupLoop::endIter(LoopContext& ctx)
{
  std::cout << "setup end\n";
  if (!ctx.app().halfTransient())
    ctx.problem().outputStep(EXEC_FINAL);
  done();
}

////////////////// MeshRefinementLoop ////////////////////
MeshRefinementLoop::MeshRefinementLoop(const InputParameters& params, LoopContext& ctx) { }

std::string MeshRefinementLoop::name()
{
  return "mesh-refinement-loop";
}

void
MeshRefinementLoop::endIter(LoopContext& ctx)
{
  unsigned int max_steps = ctx.problem().adaptivity().getSteps();
  if (iter() < max_steps)
    ctx.problem().adaptMesh();
  else
    done();
}

////////////////// TimeLoop ////////////////////
void
validParamsTimeLoop(InputParameters& params)
{
  params.addParam<unsigned int>("num_steps",       std::numeric_limits<unsigned int>::max(),     "The number of timesteps in a transient run");
}

TimeLoop::TimeLoop(const InputParameters& params, LoopContext& ctx) :
    _num_steps(params.get<unsigned int>("num_steps")),
    _steps_taken(0),

    _t_step(ctx.problem().timeStep()),
    _time(ctx.problem().time()),
    _time_old(ctx.problem().timeOld()),
    _dt(ctx.problem().dt()),
    _dt_old(ctx.problem().dtOld()),
    _sync_times(ctx.app().getOutputWarehouse().getSyncTimes())
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

void
TimeLoop::beginIter(LoopContext& ctx)
{
  std::cout << "time begin\n";
  ctx.problem().backupMultiApps(EXEC_TIMESTEP_BEGIN);
  ctx.problem().backupMultiApps(EXEC_TIMESTEP_END);

  ctx.problem().execute(EXEC_TIMESTEP_BEGIN);
  ctx.problem().outputStep(EXEC_TIMESTEP_BEGIN);
}

void
TimeLoop::endIter(LoopContext& ctx)
{
  std::cout << "time end\n";
  if (_num_steps == 1 && ctx.failed()) {
    //ctx.console << "aborting early: solve did not converge\n";
    done();
    return;
  }
  std::cout << "spot1\n";

  ctx.problem().onTimestepEnd();
  ctx.problem().execute(EXEC_TIMESTEP_END);
  std::cout << "time end\n";

  ctx.problem().computeIndicators();
  ctx.problem().computeMarkers();
  std::cout << "spot2\n";

  ctx.problem().outputStep(EXEC_TIMESTEP_END);
  if (_num_steps >= iter())
    done();
}

////////////////// SolveLoop ////////////////////
void
validParamsSolveLoop(InputParameters& params)
{
}

SolveLoop::SolveLoop(const InputParameters& params, LoopContext& ctx)
{
}

std::string SolveLoop::name()
{
  return "solve-loop";
}

void
SolveLoop::beginIter(LoopContext& ctx)
{
  std::cout << "solve begin\n";
  ctx.unfail();
  ctx.solve();
  if (!ctx.problem().converged())
    ctx.fail("solve failed to converge");
}

void
SolveLoop::endIter(LoopContext& ctx)
{
  std::cout << "solve end\n";
  done();
}

////////////////// Xfemloop ////////////////////
void
validParamsXfemLoop(InputParameters& params)
{
  params.addParam<unsigned int>("xfem_max_iter", std::numeric_limits<unsigned int>::max(), "Max number of iterations for xfem updates.");
}

XfemLoop::XfemLoop(const InputParameters& params, LoopContext& ctx)
  : _max_update(params.get<unsigned int>("xfem_max_iter"))
{
}

std::string XfemLoop::name()
{
  return "xfem-loop";
}

void
XfemLoop::beginIter(LoopContext& ctx)
{
  // TODO: does this need to call FEProblem::timestepSetup?
}

void
XfemLoop::endIter(LoopContext& ctx)
{
  if (ctx.failed() || !ctx.problem().updateMeshXFEM() || (iter() >= _max_update))
    done();
  //else
  //  _console << "XFEM modifying mesh, repeating step"<<std::endl;
}

////////////////// PicardLoop ////////////////////
void
validParamsPicardLoop(InputParameters& params)
{
  params.addParam<unsigned int>("picard_max_iter",       std::numeric_limits<unsigned int>::max(),     "Max number of picard iterations.");
  params.addParam<Real>("picard_rel_tol", 1e-8, "The relative nonlinear residual drop to shoot for during Picard iterations.  This check is performed based on the Master app's nonlinear residual.");
  params.addParam<Real>("picard_abs_tol", 1e-50, "The absolute nonlinear residual to shoot for during Picard iterations.  This check is performed based on the Master app's nonlinear residual.");
}

PicardLoop::PicardLoop(const InputParameters& params, LoopContext& ctx)
  : _max_its(params.get<unsigned int>("picard_max_iter")),
    _abs_tol(params.get<Real>("picard_abs_tol")),
    _rel_tol(params.get<Real>("picard_rel_tol"))
{
}

std::string PicardLoop::name()
{
  return "picard-loop";
}

void
PicardLoop::beginIter(LoopContext& ctx)
{
    //_console << "\nBeginning Picard Iteration " << iter() << "\n" << std::endl;
    if (iter() == 0) // First Picard iteration - need to save off the initial nonlinear residual
    {
      _initial_norm = ctx.problem().computeResidualL2Norm();
      //_console << "Initial Picard Norm: " << _initial_norm << '\n';
    }
    else
    {
      ctx.problem().restoreMultiApps(EXEC_TIMESTEP_BEGIN);
      ctx.problem().restoreMultiApps(EXEC_TIMESTEP_END);
    }

    _begin_norm = ctx.problem().computeResidualL2Norm();
    //_console << "Picard Norm after TIMESTEP_BEGIN MultiApps: " << _begin_norm << '\n';
}

void
PicardLoop::endIter(LoopContext& ctx)
{
  if (ctx.failed())
  {
    done();
    return;
  }
  if (iter() >= _max_its)
    done();

  Real end_norm = ctx.problem().computeResidualL2Norm();
  //_console << "Picard Norm after TIMESTEP_END MultiApps: " << end_norm << '\n';
  Real max_norm = std::max(_begin_norm, end_norm);
  Real max_relative_drop = max_norm / _initial_norm;
  if (max_norm < _abs_tol || max_relative_drop < _rel_tol)
  {
    //_console << "Picard converged!" << std::endl;
    done();
  }
}
