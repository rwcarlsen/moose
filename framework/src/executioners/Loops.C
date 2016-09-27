
#include "Loops.h"

template<>
InputParameters validParams<Loops>()
{
  InputParameters params = validParams<MooseObject>();
  params.registerBase("Executioner");
  return params;
}

Loops::Loops(const InputParameters & parameters) :
    MooseObject(parameters),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(parameters, "Executioners"),
    _problem(nullptr)
{
}

void Loops::initialize(QueenOfHearts& queen) {
#ifdef LIBMESH_ENABLE_AMR
  LOOP_FROM_METHODS(queen, meshAdaptivity);
#endif
  LOOP_FROM_METHODS(queen, step);
}

bool Loops::meshAdaptivityBegin(IterInfo info)
{
  return false;
}

bool Loops::meshAdaptivityEnd(IterInfo info)
{
  return false;
}


bool Loops::stepBegin(IterInfo info)
{
  return false;
}

bool Loops::stepEnd(IterInfo info)
{
  return true;
}

