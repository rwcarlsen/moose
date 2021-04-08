#pragma once

#include "graph.h"

#include "FEProblemBase.h"

class MeshLocation
{
public:
  dag::LoopType type;
  Elem * elem = nullptr;
  FaceInfo * face = nullptr;
  Node * node = nullptr;
  unsigned int side = 0;
  BoundaryID boundary = 0;
};

using UniversalRange = StoredRange<std::vector<MeshLocation *>::const_iterator, const MeshLocation *>;

void runLoop(FEProblemBase & fep,
             std::vector<std::vector<dag::Node *>> & nodes,
             const UniversalRange & range);

/**
 * This loops over a set of mesh faces (i.e. FaceInfo objects).  Callback
 * routines are provided for visiting each face, for visiting boundary faces,
 * for sudomain changes, and pre/post many of these events.
 */
class UniversalLoop
{
public:
  UniversalLoop(FEProblemBase & fe_problem, std::vector<dag::Node *> & nodes);

  UniversalLoop(UniversalLoop & x, Threads::split /*split*/)
    : _fe_problem(x._fe_problem), _mesh(x._mesh), _objects(x._objects){};

  virtual ~UniversalLoop(){};

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

