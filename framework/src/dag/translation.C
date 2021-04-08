#include "translation.h"

#include "FEProblemBase.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"
#include "MooseMesh.h"
#include "KernelBase.h"
#include "NodalBCBase.h"
#include "IntegratedBCBase.h"
#include "TimeIntegrator.h"
#include "SwapBackSentinel.h"

namespace translation
{

void
buildMeshLocations(MooseMesh & mesh,
                   std::vector<std::unique_ptr<MeshLocation>> & all,
                   std::map<dag::LoopType, std::vector<MeshLocation *>> & locs)
{
  auto begin = mesh.getMesh().active_elements_begin();
  auto end = mesh.getMesh().active_elements_end();

  for (auto it = begin; it != end; ++it)
  {
    Elem * elem = *it;
    auto block = elem->subdomain_id();
    all.emplace_back(new MeshLocation{
        dag::LoopType(dag::LoopCategory::Elemental_onElem, block), elem, nullptr, nullptr, 0, 0});

    auto loc = all.back().get();
    locs[loc->type].push_back(loc);

    for (unsigned int side = 0; side < elem->n_sides(); side++)
    {
      for (auto boundary : mesh.getBoundaryIDs(elem, side))
      {
        all.emplace_back(
            new MeshLocation{dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary),
                             elem,
                             nullptr,
                             nullptr,
                             side,
                             boundary});
        auto loc = all.back().get();
        locs[loc->type].push_back(loc);
      }
    }
  }

  ConstBndNodeRange & bnd_nodes = *mesh.getBoundaryNodeRange();
  for (const auto & bnode : bnd_nodes)
  {
    BoundaryID boundary = bnode->_bnd_id;
    Node * node = bnode->_node;
    all.emplace_back(new MeshLocation{dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary),
                                      nullptr,
                                      nullptr,
                                      node,
                                      0,
                                      boundary});
    auto loc = all.back().get();
    locs[loc->type].push_back(loc);
  }
}

void
buildSpecialNodes(FEProblemBase & fe, GraphData & gd, const std::set<TagID> & tags)
{
  for (auto block : fe.mesh().meshSubdomains())
  {
    auto elem_setup = gd.graph.create(
        "elem_setup", false, false, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
    gd.elem_setup[block] = elem_setup;
    elem_setup->setRunFunc([&fe](const MeshLocation & loc, THREAD_ID tid) {
      fe.prepare(loc.elem, tid);
      fe.reinitElem(loc.elem, tid);
    });
    auto elem_teardown = gd.graph.create(
        "elem_teardown", false, false, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
    gd.elem_teardown[block] = elem_teardown;
    elem_teardown->setRunFunc([&fe](const MeshLocation &, THREAD_ID tid) {
      fe.cacheResidual(tid);
      {
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
        fe.addCachedResidual(tid);
      }
    });
  }

  for (auto boundary : fe.mesh().meshBoundaryIds())
  {
    auto side_setup =
        gd.graph.create("side_setup",
                        false,
                        false,
                        dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
    gd.side_setup[boundary] = side_setup;
    side_setup->setRunFunc([&fe](const MeshLocation & loc, THREAD_ID tid) {
      fe.prepare(loc.elem, tid);
      fe.reinitElem(loc.elem, tid);
      fe.reinitElemFace(loc.elem, loc.side, loc.boundary, tid);
    });

    auto side_teardown =
        gd.graph.create("side_teardown",
                        false,
                        false,
                        dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
    gd.side_teardown[boundary] = side_teardown;
    side_teardown->setRunFunc([&fe](const MeshLocation &, THREAD_ID tid) {
      fe.cacheResidual(tid);
      {
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
        fe.addCachedResidual(tid);
      }
    });
    auto side_node_setup =
        gd.graph.create("side_node_setup",
                        false,
                        false,
                        dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
    gd.side_node_setup[boundary] = side_node_setup;
    side_node_setup->setRunFunc([&fe](const MeshLocation & loc, THREAD_ID tid) {
      fe.reinitNodeFace(loc.node, loc.boundary, tid);
    });
  }

  gd.residual_setup =
      gd.graph.create("residual_setup", true, false, dag::LoopType(dag::LoopCategory::None));
  gd.residual_setup->setRunFunc([&fe, tags](const MeshLocation &, THREAD_ID) {
    fe.setCurrentlyComputingResidual(true);
    fe.setCurrentExecuteOnFlag(EXEC_LINEAR);
    fe.getNonlinearSystemBase().zeroTaggedVectors(tags);
    fe.getNonlinearSystemBase().residualSetup();
  });

  gd.residual_teardown =
      gd.graph.create("residual_teardown", true, false, dag::LoopType(dag::LoopCategory::None));
  gd.residual_teardown->setRunFunc([&fe, tags](const MeshLocation &, THREAD_ID) {
    if (fe.getNonlinearSystemBase().haveResidualTimeVector())
      fe.getNonlinearSystemBase().getResidualTimeVector().close();
    fe.getNonlinearSystemBase().getResidualNonTimeVector().close();
    fe.getNonlinearSystemBase().closeTaggedVectors(tags);
    fe.setCurrentExecuteOnFlag(EXEC_NONE);
    fe.getNonlinearSystemBase().activeAllMatrixTags();
    fe.setCurrentlyComputingResidual(false);
  });

  gd.solution = gd.graph.create("solution", true, false, dag::LoopType(dag::LoopCategory::None));
  gd.solution->needs(gd.residual_teardown);

  gd.pre_nodal_residual =
      gd.graph.create("pre_nodal_residual", true, false, dag::LoopType(dag::LoopCategory::None));
  gd.pre_nodal_residual->setRunFunc([&fe, tags](const MeshLocation &, THREAD_ID) {
    fe.getNonlinearSystemBase().closeTaggedVectors(tags);
    bool required_residual =
        tags.find(fe.getNonlinearSystemBase().residualVectorTag()) == tags.end() ? false : true;
    if (required_residual)
    {
      auto & residual =
          fe.getNonlinearSystemBase().getVector(fe.getNonlinearSystemBase().residualVectorTag());
      if (fe.getNonlinearSystemBase().getTimeIntegrator())
        fe.getNonlinearSystemBase().getTimeIntegrator()->postResidual(residual);
      else
        residual += fe.getNonlinearSystemBase().getResidualNonTimeVector();
      residual.close();
    }
  });

  for (auto block : fe.mesh().meshSubdomains())
    gd.pre_nodal_residual->needs(gd.elem_teardown[block]);
}

dag::Node *
convertKernel(FEProblemBase & fe,
              GraphData & gd,
              std::vector<KernelBase *> & kernels,
              SubdomainID block)
{
  auto objname = kernels[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
  obj->setRunFunc([kernels, &fe](const MeshLocation & loc, THREAD_ID tid) {
    auto props = kernels[tid]->getMatPropDependencies();
    auto & vars = kernels[tid]->getMooseVariableDependencies();
    fe.setActiveElementalMooseVariables(vars, tid);
    fe.setActiveMaterialProperties(props, tid);
    fe.prepareMaterials(loc.type.block, tid);

    SwapBackSentinel sentinel(fe, &FEProblem::swapBackMaterials, tid);
    fe.reinitMaterials(loc.type.block, tid);

    kernels[tid]->computeResidual();
  });

  obj->needs(gd.elem_setup[block]);
  gd.elem_teardown[block]->needs(obj);
  obj->needs(gd.residual_setup);
  gd.residual_teardown->needs(obj);
  return obj;
}

dag::Node *
convertBC(FEProblemBase & fe,
          GraphData & gd,
          std::vector<IntegratedBCBase *> & bcs,
          BoundaryID boundary)
{
  auto objname = bcs[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
  obj->setRunFunc([bcs, &fe](const MeshLocation & loc, THREAD_ID tid) {
    auto props = bcs[tid]->getMatPropDependencies();
    auto & vars = bcs[tid]->getMooseVariableDependencies();
    fe.setActiveElementalMooseVariables(vars, tid);
    fe.setActiveMaterialProperties(props, tid);
    fe.prepareMaterials(loc.elem->subdomain_id(), tid);

    SwapBackSentinel sentinel(fe, &FEProblem::swapBackMaterialsFace, tid);
    fe.reinitMaterialsFace(loc.elem->subdomain_id(), tid);
    fe.reinitMaterialsBoundary(loc.boundary, tid);

    if (bcs[tid]->shouldApply())
      bcs[tid]->computeResidual();
  });

  obj->needs(gd.side_setup[boundary]);
  gd.side_teardown[boundary]->needs(obj);
  obj->needs(gd.residual_setup);
  gd.residual_teardown->needs(obj);
  return obj;
}

dag::Node *
convertNodalBC(FEProblemBase &,
               GraphData & gd,
               std::vector<NodalBCBase *> & bcs,
               BoundaryID boundary)
{
  auto objname = bcs[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
  obj->setRunFunc([bcs](const MeshLocation & /*loc*/, THREAD_ID tid) {
    if (bcs[tid]->shouldApply())
      bcs[tid]->computeResidual();
  });

  gd.residual_nodes.push_back(obj);
  obj->needs(gd.side_node_setup[boundary]);
  obj->needs(gd.pre_nodal_residual);
  gd.residual_teardown->needs(obj);
  return obj;
}

void
buildLoops(FEProblemBase & fe, const std::set<TagID> & tags, GraphData & gd)
{
  translation::buildSpecialNodes(fe, gd, tags);

  const unsigned int tmp_tid = 0;
  const auto n_threads = libMesh::n_threads();
  auto & kern_warehouse =
      fe.getNonlinearSystemBase().getKernelWarehouse().getVectorTagsObjectWarehouse(tags, tmp_tid);
  auto & ibc_warehouse =
      fe.getNonlinearSystemBase().getIntegratedBCWarehouse().getVectorTagsObjectWarehouse(tags,
                                                                                          tmp_tid);
  auto & nbc_warehouse =
      fe.getNonlinearSystemBase().getNodalBCWarehouse().getVectorTagsObjectWarehouse(tags, tmp_tid);

  for (auto block : fe.mesh().meshSubdomains())
  {
    if (!kern_warehouse.hasActiveBlockObjects(block, tmp_tid))
      continue;

    auto & kernels = kern_warehouse.getBlockObjects(block, tmp_tid);
    for (auto kern : kernels)
    {
      auto k = kern.get();
      std::vector<KernelBase *> kernel_copies;
      for (unsigned int i = 0; i < n_threads; i++)
        kernel_copies.push_back(kern_warehouse.getObject(k->name(), i).get());
      convertKernel(fe, gd, kernel_copies, block);
    }
  }
  for (auto boundary : fe.mesh().meshBoundaryIds())
  {
    if (ibc_warehouse.hasActiveBoundaryObjects(boundary, tmp_tid))
    {
      auto & ibcs = ibc_warehouse.getBoundaryObjects(boundary, tmp_tid);
      for (auto bc : ibcs)
      {
        auto b = bc.get();
        std::vector<IntegratedBCBase *> ibc_copies;
        for (unsigned int i = 0; i < n_threads; i++)
          ibc_copies.push_back(ibc_warehouse.getObject(b->name(), i).get());
        convertBC(fe, gd, ibc_copies, boundary);
      }
    }

    if (nbc_warehouse.hasActiveBoundaryObjects(boundary, tmp_tid))
    {
      auto & nbcs = nbc_warehouse.getBoundaryObjects(boundary, tmp_tid);
      for (auto bc : nbcs)
      {
        auto b = bc.get();
        std::vector<NodalBCBase *> nbc_copies;
        for (unsigned int i = 0; i < n_threads; i++)
          nbc_copies.push_back(nbc_warehouse.getObject(b->name(), i).get());
        convertNodalBC(fe, gd, nbc_copies, boundary);
      }
    }
  }

  // calculate the loops
  std::set<dag::Node *> start_nodes = {gd.solution};
  auto partitions = dag::computePartitions(gd.graph, true);
  for (auto & g : partitions)
    if (g.reachable(start_nodes))
      gd.partitions.push_back(g);
    else
      for (auto s : start_nodes)
        if (g.contains(s))
        {
          gd.partitions.push_back(g);
          break;
        }

  gd.objs = computeLoops(gd.partitions);
  for (auto nodes : gd.objs)
    gd.loop_type.push_back(nodes[0][0]->loopType());
  dag::printLoops(gd.objs);

  // show all the given subgraphs on a single graph.
  dag::Subgraph full;
  for (auto & p : gd.partitions)
    for (auto n : p.nodes())
      full.add(n);
  std::cout << dag::dotGraph(full);
}

} // namespace translation
