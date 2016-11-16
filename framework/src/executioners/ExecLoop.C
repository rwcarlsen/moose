
#include "ExecLoop.h"

LoopContext::LoopContext(MooseApp& app, FEProblem& prob)
{
  
}

MooseApp&
LoopContext::app()
{
  return _app;
}

FEProblem&
LoopContext::problem()
{
  return _prob;
}

void
LoopContext::fail(std::string reason)
{
  _failed = true;
  _failed_reason = reason;
}

void
LoopContext::unfail()
{
  _failed = false;
  _failed_reason.clear();
}

bool
LoopContext::failed()
{
  return _failed;
}

Real
LoopContext::solnDiffNorm()
{
  return _soln_diff_norm;
}

LoopContext::solve()
{
  _problem.timestepSetup();
  _problem.updateActiveObjects();

  timeval solve_start;
  timeval solve_end;
  gettimeofday(&solve_start, nullptr);

  _prob.solve();
  
  gettimeofday(&solve_end, nullptr);
  Real _solve_time = (static_cast<Real>(solve_end.tv_sec  - solve_start.tv_sec) +
                                             static_cast<Real>(solve_end.tv_usec - solve_start.tv_usec)*1.e-6);
  if (!_prob.converged())
    fail("solve failed to converge");
  _soln_diff_norm = _problem.relativeSolutionDifferenceNorm();
}

ExecLoop::~ExecLoop() { };

void
ExecLoop::addChild(ExecLoop* loop)
{
  _children.push_back(loop);
}

void
ExecLoop::run(LoopContext* ctx)
{
  runLoop(ctx, 0);
}

int
ExecLoop::iter()
{
  return iter(_iters.size() - 1);
}

int
ExecLoop::iter(int loop)
{
  if (loop >= _iters.size())
    return 0;
  return _iters[loop];
}

int
ExecLoop::iter(std::string loop)
{
  for (int i = 0; i < _iters.size(); i++)
  {
    if (_names[i] == loop)
      return _iters[i];
  }
  mooseError("iteration number requested for invalid loop name '" + loop + "'");
}

void
ExecLoop::runLoop(LoopContext* ctx, int loop)
{
  _iters.push_back(0);
  _names.push_back(name());
  while (true)
  {
    _iters[loop]++;
    
    _done = false;
    beginIter(ctx);
    if (_done)
      break;

    for (auto child : _children)
    {
      child->_iters = _iters;
      child->runLoop(ctx, loop + 1);
    }
    
    _done = false;
    endIter(ctx);
    if (_done)
      break;
  }
  _iters.pop_back();
  _names.pop_back();

  // reset iter counts
  for (int i = loop; i < _iters.size(); i++) _iters[i] = 0;
}

void
LoopContext::done()
{
  _done = true;
}

