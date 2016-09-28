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

class FEProblem;
class MooseApp;
class Console;

struct LoopContext {
  MooseApp* app;
  FEProblem* problem;
};

class ExecLoop {
public:
  virtual ~ExecLoop();

  void run(LoopContext& ctx);
  void addChild(ExecLoop* loop);

  int iter();
  int iter(int loop);
  int iter(std::string loop);

  virtual std::string name() = 0;
  virtual bool beginIter(LoopContext& ctx) = 0;
  virtual bool endIter(LoopContext& ctx) = 0;

private:
  void runLoop(LoopContext& ctx, int loop);

  std::vector<int> _iters;
  std::vector<std::string> _names;
  std::vector<ExecLoop*> _children;
};

#endif //EXECLOOP_H
