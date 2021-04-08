
#include "loop.h"

void
runLoop(FEProblemBase & fep,
        std::vector<std::vector<dag::Node *>> & nodes,
        const UniversalRange & range)
{
  std::vector<dag::Node *> all;
  for (auto group : nodes)
    for (auto n : group)
      all.push_back(n);

  UniversalLoop loop(fep, all);
  Threads::parallel_reduce(range, loop);
  for (auto n : all)
    n->join();
}

UniversalLoop::UniversalLoop(FEProblemBase & fe_problem, std::vector<dag::Node *> & nodes)
  : _fe_problem(fe_problem), _mesh(fe_problem.mesh()), _objects(nodes){};

// TODO: ensure the vector<faceinfo> data structure needs to be built such
// that for all sides on an interface between two subdomains, the elements of
// the same subdomain are used consistently for all the "elem" (i.e. not
// "neighbor") parameters in order to avoid jumping back and forth along the
// boundary between using one or the other subdomains' FV kernels
// unpredictably.
void
UniversalLoop::operator()(const UniversalRange & range, bool bypass_threading)
{
  try
  {
    try
    {
      ParallelUniqueId puid;
      _tid = bypass_threading ? 0 : puid.id;

      UniversalRange::const_iterator loc = range.begin();
      for (loc = range.begin(); loc != range.end(); ++loc)
        for (auto obj : _objects)
          obj->run(**loc, _tid);
    }
    catch (libMesh::LogicError & e)
    {
      throw MooseException("We caught a libMesh error");
    }
  }
  catch (MooseException & e)
  {
    caughtMooseException(e);
  }
}
