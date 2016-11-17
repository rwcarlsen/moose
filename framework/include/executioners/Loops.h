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
class StepperBlock;

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
  void setupTimeIntegrator(const InputParameters& params, LoopContext& ctx);

  std::unique_ptr<StepperBlock> _stepper;
  unsigned int _num_steps;
  MooseEnum _time_scheme;
  Real _time;
  Real _start_time;
  Real _end_time;
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
  Real _abs_tol;
  Real _rel_tol;
  Real _initial_norm;
  Real _begin_norm;
};

#endif //LOOPS_H

