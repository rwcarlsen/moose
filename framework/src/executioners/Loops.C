
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
  LOOP_FROM_METHODS(queen, MeshAdaptivity);
#endif
  LOOP_FROM_METHODS(queen, Step);
}

bool Loops::beginMeshAdaptivity(IterInfo info)
{
}

bool Loops::endMeshAdaptivity(IterInfo info)
{
}


bool Loops::beginStep(IterInfo info)
{
}

bool Loops::endStep(IterInfo info)
{
}

