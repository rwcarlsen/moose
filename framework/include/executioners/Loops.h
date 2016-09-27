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

  void initialize(QueenOfHearts& queen);

  bool meshAdaptivityBegin(IterInfo info);
  bool meshAdaptivityEnd(IterInfo info);

  bool stepBegin(IterInfo info);
  bool stepEnd(IterInfo info);

private:
  FEProblem* _problem;

};

#endif //LOOPS_H

