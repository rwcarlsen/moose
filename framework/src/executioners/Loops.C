
#include "Loops.h"
#include "Stepper.h"
#include "TimeStepper.h"
#include "Predictor.h"

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
  
  _root.reset(new SetupLoop(_pars, ctx));
  ExecLoop* curr = _root.get();

  std::string flavor = _pars.get<std::string>("flavor");
  if (flavor == "steady-state")
  {
#ifdef LIBMESH_ENABLE_AMR
    curr = curr->addChild(new MeshRefinementLoop(_pars, ctx));
#endif
    curr = curr->addChild(new TimeLoop(_pars, ctx));
    curr = curr->addChild(new SolveLoop(_pars, ctx));
  }
  else if (flavor == "transient")
  {
    curr = curr->addChild(new TimeLoop(_pars, ctx));
    curr = curr->addChild(new SolveLoop(_pars, ctx));
  }
  else
  {
    mooseError("unrecognized Loops flavor '" + flavor + "'");
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
MeshRefinementLoop::MeshRefinementLoop(const InputParameters& params, LoopContext& ctx)
{
  _max_steps = ctx.problem().adaptivity().getSteps();
}

std::string MeshRefinementLoop::name()
{
  return "mesh-refinement-loop";
}

void
MeshRefinementLoop::endIter(LoopContext& ctx)
{
  if (ctx.problem().adaptivity().isOn())
    ctx.problem().adaptMesh();
  if (iter() >= _max_steps)
    done();
}

////////////////// TimeLoop ////////////////////
void
validParamsTimeLoop(InputParameters& params)
{
  params.addParam<unsigned int>("num_steps",       std::numeric_limits<unsigned int>::max(),     "The number of timesteps in a transient run");
  MooseEnum schemes("implicit-euler explicit-euler crank-nicolson bdf2 rk-2 dirk explicit-tvd-rk-2", "implicit-euler");
  params.addParam<MooseEnum>("scheme", schemes, "Time integration scheme used.");
  params.addParam<Real>("start_time",      0.0,    "The start time of the simulation");
  params.addParam<Real>("end_time",        1.0e30, "The end time of the simulation");
  params.addParam<bool>("abort_on_solve_fail", false, "abort if solve not converged rather than cut timestep");
  params.addParam<Real>("timestep_tolerance", 2.0e-14, "the tolerance setting for final timestep size and sync times");
  params.addParam<Real>("dt",              1.,     "The timestep size between solves");
  params.addParam<Real>("dtmin",           2.0e-14,    "The minimum timestep size in an adaptive run");
  params.addParam<Real>("dtmax",           1.0e30, "The maximum timestep size in an adaptive run");
}

TimeLoop::TimeLoop(const InputParameters& params, LoopContext& ctx)
  : _stepper(nullptr),
    _si(),
    _num_steps(params.get<unsigned int>("num_steps")),
    _time_scheme(params.get<MooseEnum>("scheme")),
    _tol(params.get<Real>("timestep_tolerance")),
    _time(params.get<Real>("start_time")),
    _t_step(0),
    _start_time(params.get<Real>("start_time")),
    _end_time(params.get<Real>("end_time")),
    _dtmin(params.get<Real>("dtmin")),
    _dtmax(params.get<Real>("dtmax")),
    _prev_dt(0),
    _abort(params.get<bool>("abort_on_solve_fail")),
    _last_converged(true),
    _snapshots()
{
  std::string flavor = params.get<std::string>("flavor");
  if (flavor == "steady-state")
  {
    _num_steps = 1;
    return;
  }
  
  setupTimeStepper(params, ctx);
  setupTimeIntegrator(params, ctx);

  ctx.problem().transient(true);
  // Either a start_time has been forced on us, or we want to tell the App about what our start time is (in case anyone else is interested.
  if (ctx.app().hasStartTime())
  {
    _start_time = ctx.app().getStartTime();
    _time = _start_time;
  }
  else if (params.isParamSetByUser("start_time"))
    ctx.app().setStartTime(_start_time);

  if (ctx.app().halfTransient()) // Cut timesteps and end_time in half...
  {
    _end_time /= 2.0;
    _num_steps /= 2.0;
    if (_num_steps == 0) // Always do one step in the first half
      _num_steps = 1;
  }
}

std::string TimeLoop::name()
{
  return "time-loop";
}

void
TimeLoop::beginIter(LoopContext& ctx)
{
  std::cout << "time begin\n";
  if (!ctx.failed())
  {
    _t_step++;
#ifdef LIBMESH_ENABLE_AMR
    if (ctx.problem().adaptivity().isOn())
      ctx.problem().adaptMesh();
#endif
  }
  ctx.unfail();
  ctx.problem().timeStep() = _t_step;

  ctx.problem().backupMultiApps(EXEC_TIMESTEP_BEGIN);
  ctx.problem().backupMultiApps(EXEC_TIMESTEP_END);

  ctx.problem().onTimestepBegin();

  ctx.problem().dt() = calcDT(ctx);
  ctx.problem().time() = _time + ctx.problem().dt();
  
  bool auto_advance = false; // always done in endIter
  ctx.problem().execTransfers(EXEC_TIMESTEP_BEGIN);
  bool multi_converged = ctx.problem().execMultiApps(EXEC_TIMESTEP_BEGIN, auto_advance);
  if (!multi_converged)
  {
    ctx.fail("multi apps failed to converge");
    return;
  }

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
  else if (iter() >= _num_steps || (std::abs(_time - _end_time) <= _tol))
    done();
  else if (ctx.failed() && _abort)
    done();
  
  std::cout << "spot1\n";

  retrieveSolveState(ctx);
  std::cout << "converged = " << _last_converged << "\n";
  if (ctx.failed())
  {
    ctx.problem().restoreMultiApps(EXEC_TIMESTEP_BEGIN, true);
    ctx.problem().restoreMultiApps(EXEC_TIMESTEP_END, true);
  }
  else
  {
    ctx.problem().advanceMultiApps(EXEC_TIMESTEP_BEGIN);
    ctx.problem().advanceMultiApps(EXEC_TIMESTEP_END);
    _time += ctx.problem().dt();
  }
  
  ctx.problem().onTimestepEnd();
  ctx.problem().execute(EXEC_TIMESTEP_END);
  std::cout << "spot2\n";

  ctx.problem().computeIndicators();
  ctx.problem().computeMarkers();
  std::cout << "spot3\n";

  ctx.problem().outputStep(EXEC_TIMESTEP_END);

  if (!ctx.failed())
    ctx.problem().advanceState();
}

void
TimeLoop::setupTimeStepper(const InputParameters& params, LoopContext& ctx)
{
  StepperBlock * inner = nullptr;
  if (ctx.app()._time_stepper)
   inner = ctx.app()._time_stepper->buildStepper();
  else
  {
    if (!params.isParamSetByAddParam("end_time") && 
        !params.isParamSetByAddParam("num_steps") && params.isParamSetByAddParam("dt"))
      inner = BaseStepper::constant((_end_time - _time) / _num_steps);
    else
      inner = BaseStepper::constant(params.get<Real>("dt"));
  }
  
  std::vector<Real> sync_times;
  for (auto val : ctx.app().getOutputWarehouse().getSyncTimes())
    sync_times.push_back(val);

  // these are global/sim constraints for *every* time stepper:
  inner = BaseStepper::dtLimit(inner, _dtmin, _dtmax);
  if (sync_times.size() > 0)
    inner = BaseStepper::min(BaseStepper::fixedTimes(sync_times, _tol), inner, _tol);
  if (!ctx.app().halfTransient())
    inner = BaseStepper::bounds(inner, _start_time, _end_time);
  _stepper.reset(inner);
}

void
TimeLoop::setupTimeIntegrator(const InputParameters& params, LoopContext& ctx)
{
  if (ctx.problem().hasTimeIntegrator() && params.isParamSetByUser("scheme"))
    mooseError("You cannot specify time_scheme in the Executioner and independently add a TimeIntegrator to the system at the same time");
  if (!_time_scheme.isValid())
    mooseError("Unknown scheme");

  std::string ti_str;
  switch (_time_scheme)
  {
  case 0: ti_str = "ImplicitEuler"; break;
  case 1: ti_str = "ExplicitEuler"; break;
  case 2: ti_str = "CrankNicolson"; break;
  case 3: ti_str = "BDF2"; break;
  case 4: ti_str = "ExplicitMidpoint"; break;
  case 5: ti_str = "LStableDirk2"; break;
  case 6: ti_str = "ExplicitTVDRK2"; break;
  }

  InputParameters pars = ctx.app().getFactory().getValidParams(ti_str);
  ctx.problem().addTimeIntegrator(ti_str, ti_str, pars);
}

Real
TimeLoop::calcDT(LoopContext& ctx)
{
  updateStepperInfo(ctx);
  Real dt = _stepper->next(_si);
  if (_si.wantSnapshot())
    // the time used in this key must be *exactly* that on the stepper just saw in StepperInfo
    _snapshots[_si.time()] = ctx.app().backup();
  if (_si.rewindTime() != -1)
  {
    if (_snapshots.count(_si.rewindTime()) == 0)
      mooseError("no snapshot available for requested rewind time");
    ctx.app().restore(_snapshots[_si.rewindTime()]);
    // second calls necessary because rewind modifies state _si depends on
    updateStepperInfo(ctx);
    dt = _stepper->next(_si);
  }
  return dt;
}

void
TimeLoop::retrieveSolveState(LoopContext& ctx)
{
  _soln_nonlin.resize(ctx.problem().getNonlinearSystem().currentSolution()->size());
  _soln_aux.resize(ctx.problem().getAuxiliarySystem().currentSolution()->size());
  _soln_predicted.resize(ctx.problem().getNonlinearSystem().currentSolution()->size());
  ctx.problem().getNonlinearSystem().currentSolution()->localize(_soln_nonlin);
  ctx.problem().getAuxiliarySystem().currentSolution()->localize(_soln_aux);
  Predictor* p = ctx.problem().getNonlinearSystem().getPredictor();
  if (p)
  {
    _soln_predicted.resize(p->solutionPredictor().size());
    p->solutionPredictor().localize(_soln_predicted);
  }

  _nl_its = ctx.problem().getNonlinearSystem().nNonlinearIterations();
  _l_its = ctx.problem().getNonlinearSystem().nLinearIterations();
  _last_converged = ctx.problem().converged();
}

void
TimeLoop::updateStepperInfo(LoopContext& ctx)
{
  _soln_nonlin.resize(ctx.problem().getNonlinearSystem().currentSolution()->size());
  _soln_aux.resize(ctx.problem().getAuxiliarySystem().currentSolution()->size());
  _soln_predicted.resize(ctx.problem().getNonlinearSystem().currentSolution()->size());

  if (_si.dt() != _prev_dt || iter() == 0)
    _si.pushHistory(_prev_dt, true, 0);

  _si.update(iter(), _time, ctx.problem().dt(), _nl_its, _l_its, _last_converged,
             ctx.solveTime(), _soln_nonlin, _soln_aux, _soln_predicted);

  _prev_dt = _si.dt(); // for restart
};

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
    _rel_tol(params.get<Real>("picard_rel_tol")),
    _initial_norm(0),
    _begin_norm(0)
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
