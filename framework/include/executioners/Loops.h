#ifndef LOOPS_H
#define LOOPS_H

#include "InputParameters.h"
#include "QueenOfHearts.h"

// System includes
#include <string>

class Loops

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

  bool beginMeshAdaptivity(IterInfo info);
  bool endMeshAdaptivity(IterInfo info);

  bool beginStep(IterInfo info);
  bool endStep(IterInfo info);

private:
  FEProblem* _problem;

};

#endif //LOOPS_H

