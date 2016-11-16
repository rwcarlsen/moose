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

#ifndef EXECLOOP_H
#define EXECLOOP_H

#include <vector>
#include <string>

#include "Moose.h"

class FEProblem;
class MooseApp;
class Console;

class LoopContext {
public:
  LoopContext(MooseApp& app, FEProblem& prob);
  MooseApp& app();
  FEProblem& problem();

  /// Returns true if the problem solve converged.
  void solve();
  void fail(std::string reason = "none");
  void unfail();
  bool failed();
  std::string failedReason();
  Real solveTime();

private:
  MooseApp& _app;
  FEProblem& _prob;
  Real _solve_time;
  bool _failed;
  std::string _failed_reason;
};

class ExecLoop {
public:
  virtual ~ExecLoop();

  void run(LoopContext& ctx);
  ExecLoop* addChild(ExecLoop* loop);
  void done();

  int iter();
  int iter(int loop);
  int iter(std::string loop);

  virtual std::string name() = 0;
  virtual void beginIter(LoopContext& ctx) {};
  virtual void endIter(LoopContext& ctx) {};

private:
  void runLoop(LoopContext& ctx, int loop);

  std::vector<int> _iters;
  std::vector<std::string> _names;
  std::vector<ExecLoop*> _children;
  bool _done;
};

#endif //EXECLOOP_H
