#pragma once

#include "loop.h"

class KernelBase;
class NodalBCBase;
class IntegratedBCBase;

namespace translation
{

struct GraphData
{
  dag::Graph graph;
  std::vector<dag::Subgraph> partitions;
  std::vector<std::vector<std::vector<dag::Node *>>> objs;
  std::map<SubdomainID, dag::Node *> elem_setup;
  std::map<SubdomainID, dag::Node *> elem_teardown;
  std::map<BoundaryID, dag::Node *> side_setup;
  std::map<BoundaryID, dag::Node *> side_teardown;
  std::map<BoundaryID, dag::Node *> side_node_setup;
  std::vector<dag::Node *> residual_nodes;
  std::vector<dag::LoopType> loop_type;
  dag::Node * residual_setup = nullptr;
  dag::Node * residual_teardown = nullptr;
  dag::Node * pre_nodal_residual = nullptr;
  dag::Node * solution = nullptr;
};

void buildMeshLocations(MooseMesh & mesh,
                        std::vector<std::unique_ptr<MeshLocation>> & all,
                        std::map<dag::LoopType, std::vector<MeshLocation *>> & locs);
void buildLoops(FEProblemBase & fe, const std::set<TagID> & tags, GraphData & gd);

dag::Node * convertKernel(FEProblemBase & fe,
                          GraphData & gd,
                          std::vector<KernelBase *> & kernels,
                          SubdomainID block);
dag::Node * convertBC(FEProblemBase & fe,
                      GraphData & gd,
                      std::vector<IntegratedBCBase *> & bcs,
                      BoundaryID boundary);
dag::Node * convertNodalBC(FEProblemBase &,
                           GraphData & gd,
                           std::vector<NodalBCBase *> & bcs,
                           BoundaryID boundary);

} // namespace translation
