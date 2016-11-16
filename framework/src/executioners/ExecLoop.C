
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

bool
LoopContext::solve()
{
  _prob.solve();
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
    if (beginIter(ctx)) break;
    for (auto child : _children)
    {
      child->_iters = _iters;
      child->runLoop(ctx, loop + 1);
    }
    if (endIter(ctx)) break;
  }
  _iters.pop_back();
  _names.pop_back();

  // reset iter counts
  for (int i = loop; i < _iters.size(); i++) _iters[i] = 0;
}

