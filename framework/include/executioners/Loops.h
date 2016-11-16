#ifndef LOOPS_H
#define LOOPS_H

#include "MooseObject.h"
#include "UserObjectInterface.h"
#include "PostprocessorInterface.h"
#include "Restartable.h"
#include "InputParameters.h"
#include "ExecLoop.h"

// System includes
#include <string>

class Loops;

template<>
InputParameters validParams<Loops>();

class Loops :
  public MooseObject,
  public UserObjectInterface,
  public PostprocessorInterface,
  public Restartable
{
public:
  Loops(const InputParameters & parameters);

  void initialize(FEProblem* problem);
  void run();

private:
  FEProblem* _problem;
  ExecLoop* _root;
};

class SetupLoop : public ExecLoop
{
public:
  SetupLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void beginIter(LoopContext& ctx);
  virtual void endIter(LoopContext& ctx);
private:
  bool _steady;
};

class MeshRefinementLoop : public ExecLoop
{
public:
  MeshRefinementLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void endIter(LoopContext& ctx);
};

class TimeLoop : public ExecLoop
{
public:
  TimeLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void beginIter(LoopContext& ctx);
  virtual void endIter(LoopContext& ctx);

private:
  unsigned int _num_steps;
  unsigned int _steps_taken;

  int & _t_step;
  /// Current time
  Real & _time;
  /// Previous time
  Real & _time_old;
  /// Current delta t... or timestep size.
  Real & _dt;
  Real & _dt_old;
  std::set<Real> & _sync_times;
};

class SolveLoop : public ExecLoop
{
public:
  SolveLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void beginIter(LoopContext& ctx);
  virtual void endIter(LoopContext& ctx);
};

class XfemLoop : public ExecLoop
{
public:
  XfemLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void beginIter(LoopContext& ctx);
  virtual void endIter(LoopContext& ctx);
  
private:
  int _max_update;
};

class PicardLoop : public ExecLoop
{
public:
  PicardLoop(const InputParameters& params, LoopContext& ctx);
  virtual std::string name();
  virtual void beginIter(LoopContext& ctx);
  virtual void endIter(LoopContext& ctx);
  
private:
  int _max_its;
  Real _begin_norm;
  Real _abs_tol;
  Real _rel_tol;
};

#endif //LOOPS_H

