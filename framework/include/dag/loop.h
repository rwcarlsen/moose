#include "graph.h"

#include "FEProblemBase.h"

class MeshLocation
{
public:
  dag::LoopType type;
  SubdomainID block;
  Elem * elem;
  FaceInfo * face;
  Node * node;
};

using UniversalRange = StoredRange<std::vector<MeshLocation *>::const_iterator, const MeshLocation *>;

/**
 * This loops over a set of mesh faces (i.e. FaceInfo objects).  Callback
 * routines are provided for visiting each face, for visiting boundary faces,
 * for sudomain changes, and pre/post many of these events.
 */
class UniversalLoop
{
public:
  UniversalLoop(FEProblemBase & fe_problem, std::vector<dag::Node *> & nodes);

  UniversalLoop(UniversalLoop & x, Threads::split split);

  virtual ~UniversalLoop();

  virtual void operator()(const UniversalRange & range, bool bypass_threading = false);

  void join(const UniversalLoop & /*y*/){};

  /// Called if a MooseException is caught anywhere during the computation.
  virtual void caughtMooseException(MooseException &){};

private:
  FEProblemBase & _fe_problem;
  MooseMesh & _mesh;
  THREAD_ID _tid;
  std::vector<dag::Node *> & _objects;
};

UniversalLoop::UniversalLoop(FEProblemBase & fe_problem,
                                              std::vector<dag::Node *> & nodes)
  : _fe_problem(fe_problem), _mesh(fe_problem.mesh()), _objects(nodes)
{
}

UniversalLoop::UniversalLoop(UniversalLoop & x, Threads::split /*split*/)
  : _fe_problem(x._fe_problem), _mesh(x._mesh), _objects(x._objects)
{
}

UniversalLoop::~UniversalLoop()
{
}

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

void runLoop(FEProblemBase & fep, std::vector<dag::Node *> & nodes, const UniversalRange & range)
{
  UniversalLoop loop(fep, nodes);
  Threads::parallel_reduce(range, loop);
  for (auto n : nodes)
    n->join();
}

