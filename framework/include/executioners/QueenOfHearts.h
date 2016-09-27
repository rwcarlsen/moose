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

#ifndef QUEENOFHEARTS_H
#define QUEENOFHEARTS_H

#include <vector>
#include <string>
#include <functional>

class FEProblem

struct IterInfo {
  std::string& name;
  const int& curr_loop_iter;
  const std::vector<int>& all_iters;
  FEProblem* problem;
};

class ExecLoop {
 public:
  virtual bool beginIter(IterInfo info) = 0;
  virtual bool endIter(IterInfo info) = 0;
};

class QueenOfHearts {
 public:
  QueenOfHearts();
  virtual ~QueenOfHearts();

  void addLoop(std::string name, ExecLoop* loop);
  void run();

  void reset();

private:
  void runLoop(int loop);

  std::vector<int> _loop_counts;
  std::vector<std::string> _loop_names;
  std::vector<ExecLoop*> _loops;
};

/////////////////////////////////////////////////////////////////////////////////
// all code below here is bonus for using labmdas and custom-named methods for //
// loop funcs all together in a single object/class                            //
/////////////////////////////////////////////////////////////////////////////////

#define LOOP_FROM_METHODS(ex,name) ex.addLoop( #name, \
        new ExecLoopFunc([this](IterInfo info){return name##Begin(info);}, \
                         [this](IterInfo info){return name##End(info);}) \
    ); \

typedef std::function< bool(IterInfo info) > LoopFunc;

class ExecLoopFunc : public ExecLoop {
 public:
  ExecLoopFunc(LoopFunc begin, LoopFunc end) : _begin(begin), _end(end) {};
  virtual bool beginIter(IterInfo info) override { return _begin(info);}
  virtual bool endIter(IterInfo info) override { return _end(info);}

 private:
  LoopFunc _begin;
  LoopFunc _end;
};

#endif //QUEENOFHEARTS_H
