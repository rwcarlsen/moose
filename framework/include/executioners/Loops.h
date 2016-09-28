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

  virtual void initialize(QueenOfHearts& queen, FEProblem* problem);

  bool setup(IterInfo info);
  bool teardown(IterInfo info);

  bool meshRefinementBegin(IterInfo info);
  bool meshRefinementEnd(IterInfo info);

  bool timeStepBegin(IterInfo info);
  bool timeStepEnd(IterInfo info);

  bool solveBegin(IterInfo info);
  bool solveEnd(IterInfo info);

private:
  FEProblem* _problem;
  std::string _flavor;
  unsigned int _max_time_steps;
};

#endif //LOOPS_H

