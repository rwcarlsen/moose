#ifndef LOOPS_H
#define LOOPS_H

#include "MooseObject.h"
#include "UserObjectInterface.h"
#include "PostprocessorInterface.h"
#include "Restartable.h"
#include "InputParameters.h"
#include "QueenOfHearts.h"

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
  std::string _flavor;
};

class SetupLoop : public ExecLoop
{
public:
  SetupLoop();
  virtual std::string name();
  virtual bool beginIter(LoopContext& ctx);
  virtual bool endIter(LoopContext& ctx);
  bool _steady;
};

class MeshRefinementLoop : public ExecLoop
{
public:
  virtual std::string name();
  virtual bool beginIter(LoopContext& ctx);
  virtual bool endIter(LoopContext& ctx);
};

class TimeLoop : public ExecLoop
{
public:
  TimeLoop();
  virtual std::string name();
  virtual bool beginIter(LoopContext& ctx);
  virtual bool endIter(LoopContext& ctx);
  bool _steady;
};

class SolveLoop : public ExecLoop
{
public:
  virtual std::string name();
  virtual bool beginIter(LoopContext& ctx);
  virtual bool endIter(LoopContext& ctx);
};

#endif //LOOPS_H

