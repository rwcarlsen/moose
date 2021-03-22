
/**
 * This loops over a set of mesh faces (i.e. FaceInfo objects).  Callback
 * routines are provided for visiting each face, for visiting boundary faces,
 * for sudomain changes, and pre/post many of these events.
 */
template <typename RangeType>
class UniversalLoop
{
public:
  UniversalLoop(FEProblemBase & fe_problem, const std::set<TagID> & tags, std::vector<DAGNode *> & nodes);

  UniversalLoop(UniversalLoop & x, Threads::split split);

  virtual ~UniversalLoop();

  virtual void operator()(const RangeType & range, bool bypass_threading = false);

  void join(const UniversalLoop & /*y*/){};

  virtual void onLocation(const MeshLocation & loc) = 0;

  /// Called if a MooseException is caught anywhere during the computation.
  virtual void caughtMooseException(MooseException &){};

private:
  FEProblemBase & _fe_problem;
  MooseMesh & _mesh;
  const std::set<TagID> & _tags;
  THREAD_ID _tid;
  std::vector<DAGNode *> & _objects;
};

template <typename RangeType>
UniversalLoop<RangeType>::UniversalLoop(FEProblemBase & fe_problem,
                                              const std::set<TagID> & tags,
                                              std::vector<DAGNode *> & nodes)
  : _fe_problem(fe_problem), _mesh(fe_problem.mesh()), _tags(tags), _objects(nodes)
{
}

template <typename RangeType>
UniversalLoop<RangeType>::UniversalLoop(ThreadedFaceLoop & x, Threads::split /*split*/)
  : _fe_problem(x._fe_problem), _mesh(x._mesh), _tags(x._tags), _objects(x._objects)
{
}

template <typename RangeType>
UniversalLoop<RangeType>::~ThreadedFaceLoop()
{
}

// TODO: ensure the vector<faceinfo> data structure needs to be built such
// that for all sides on an interface between two subdomains, the elements of
// the same subdomain are used consistently for all the "elem" (i.e. not
// "neighbor") parameters in order to avoid jumping back and forth along the
// boundary between using one or the other subdomains' FV kernels
// unpredictably.
template <typename RangeType>
void
UniversalLoop<RangeType>::operator()(const RangeType & range, bool bypass_threading)
{
  try
  {
    try
    {
      ParallelUniqueId puid;
      _tid = bypass_threading ? 0 : puid.id;

      typename RangeType::const_iterator loc = range.begin();
      for (loc = range.begin(); loc != range.end(); ++loc)
      {
        const Elem & elem = (*faceinfo)->elem();
        for (auto obj : _objects)
          obj->run(*loc, _tid);

      }
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

template <typename RangeType>
void runLoop(FEProblemBase & fep, const std::set<TagID> tags, std::vector<DAGNode *> & nodes, const RangeType & range)
{
  UniversalLoop loop(fep, tags, nodes);
  Threads::parallel_reduce(range, loop);
  for (auto n : nodes)
    n->join();
}

