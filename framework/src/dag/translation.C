#include "translation.h"

#include "FEProblemBase.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"
#include "MooseMesh.h"
#include "KernelBase.h"
#include "TimeKernel.h"
#include "NodalBCBase.h"
#include "IntegratedBCBase.h"
#include "TimeIntegrator.h"
#include "MooseVariableBase.h"
#include "Material.h"
#include "MaterialPropertyInterface.h"
#include "AuxKernel.h"

namespace translation
{

const dag::LoopType NoneType = dag::LoopType(dag::LoopCategory::None);

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

  auto & nodes = *mesh.getLocalNodeRange();
  for (const auto & node : nodes)
  {
    auto & block_ids = mesh.getNodeBlockIds(*node);
    auto block = *block_ids.begin();
    all.emplace_back(new MeshLocation{
        dag::LoopType(dag::LoopCategory::Nodal, block), nullptr, nullptr, node, 0, 0});
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
      fe.assembly(tid).setCurrentSubdomainID(loc.type.block);
      fe.prepare(loc.elem, tid);
      fe.getMaterialData(Moose::BLOCK_MATERIAL_DATA, tid)->swap(*loc.elem);
    });
    auto elem_teardown = gd.graph.create(
        "elem_teardown", false, false, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
    gd.elem_teardown[block] = elem_teardown;
    elem_teardown->setRunFunc([&fe](const MeshLocation &, THREAD_ID tid) {
      fe.swapBackMaterials(tid);
      fe.cacheResidual(tid);
      {
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
        fe.addCachedResidual(tid);
      }
    });

    auto node_setup =
        gd.graph.create("node_setup", false, false, dag::LoopType(dag::LoopCategory::Nodal, block));
    gd.node_setup[block] = node_setup;
    node_setup->setRunFunc(
        [&fe](const MeshLocation & loc, THREAD_ID tid) { fe.reinitNode(loc.node, tid); });

    // if anything ever goes in this node's run function, we'll need to add
    // some dependencies on node_teardown so it isn't pruned away.  Right now
    // nothing depends on it.
    auto node_teardown = gd.graph.create(
        "node_teardown", false, false, dag::LoopType(dag::LoopCategory::Nodal, block));
    gd.node_teardown[block] = node_teardown;
    node_teardown->setRunFunc([](const MeshLocation &, THREAD_ID) {});
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
      fe.assembly(tid).setCurrentSubdomainID(loc.elem->subdomain_id());
      fe.assembly(tid).setCurrentBoundaryID(loc.boundary);
      fe.prepare(loc.elem, tid);
      fe.reinitElemFace(loc.elem, loc.side, loc.boundary, tid);
      fe.getMaterialData(Moose::FACE_MATERIAL_DATA, tid)->swap(*loc.elem, loc.side);
    });

    auto side_teardown =
        gd.graph.create("side_teardown",
                        false,
                        false,
                        dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
    gd.side_teardown[boundary] = side_teardown;
    side_teardown->setRunFunc([&fe](const MeshLocation &, THREAD_ID tid) {
      fe.swapBackMaterialsFace(tid);
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

  gd.residual_setup = gd.graph.create("residual_setup", true, false, NoneType);
  gd.residual_setup->setRunFunc([&fe, tags](const MeshLocation &, THREAD_ID) {
    fe.getNonlinearSystemBase().zeroVariablesForResidual();
    fe.getAuxiliarySystem().zeroVariablesForResidual();

    fe.setCurrentlyComputingResidual(true);
    fe.setCurrentExecuteOnFlag(EXEC_LINEAR);
    fe.getNonlinearSystemBase().zeroTaggedVectors(tags);
    fe.getNonlinearSystemBase().residualSetup();
    fe.getAuxiliarySystem().residualSetup();
    for (unsigned int i = 0; i < libMesh::n_threads(); i++)
      fe.clearActiveElementalMooseVariables(i);
  });

  gd.residual_teardown = gd.graph.create("residual_teardown", true, false, NoneType);
  gd.residual_teardown->setRunFunc([&fe, tags](const MeshLocation &, THREAD_ID) {
    if (fe.getNonlinearSystemBase().haveResidualTimeVector())
      fe.getNonlinearSystemBase().getResidualTimeVector().close();
    fe.getNonlinearSystemBase().getResidualNonTimeVector().close();
    fe.getNonlinearSystemBase().closeTaggedVectors(tags);
    fe.setCurrentExecuteOnFlag(EXEC_NONE);
    fe.getNonlinearSystemBase().activeAllMatrixTags();
    fe.setCurrentlyComputingResidual(false);
  });

  gd.solution = gd.graph.create("solution", true, false, NoneType);
  gd.solution->needs(gd.residual_teardown);
  gd.solution->preserve();

  gd.pre_nodal_residual = gd.graph.create("pre_nodal_residual", true, false, NoneType);
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

  gd.residual_teardown->needs(gd.pre_nodal_residual);

  for (auto block : fe.mesh().meshSubdomains())
    gd.pre_nodal_residual->needs(gd.elem_teardown[block]);
  for (auto block : fe.mesh().meshBoundaryIds())
    gd.pre_nodal_residual->needs(gd.side_teardown[block]);

  gd.time_derivative = gd.graph.create("time_derivative", true, false, NoneType);
  gd.time_derivative->setRunFunc([&fe](const MeshLocation &, THREAD_ID) {
    fe.getNonlinearSystemBase().computeTimeDerivatives();
  });

  // this is not cached in order to allow it to be duplicated into every loop
  // where it is used - we want the aux finalize stuff to occur directly in
  // every loop that needs a finalized aux var otherwise we could be depending
  // on updated aux info that hasn't been finalized because the dag assumed
  // the finalization only needed to occur once which was in a loop prior to
  // the latest aux var/kernel calls/mods.
  gd.aux_soln_finalize = gd.graph.create("aux_solution_finalize", false, false, NoneType);
  gd.aux_soln_finalize->setRunFunc([&fe](const MeshLocation &, THREAD_ID) {
    auto & sys = fe.getAuxiliarySystem();
    if (!sys.isDirty())
      return;
    sys.isDirty(false);
    sys.solution().close();
    sys.system().update();
    auto ti = sys.getTimeIntegrator();
    if (fe.dt() > 0. && ti)
      ti->computeTimeDerivatives();
  });
}

// find all dag nodes corresponding to the moose variables obj depends on and
// make n (which is the dag node for the moose object obj) depend on them.
// If live vars is true, then n will depend on the "live" versions if any
// exist of the variables it depends on.
void
addVarDeps(GraphData & gd,
           MooseVariableDependencyInterface * obj,
           dag::Node * n,
           bool need_live_vars, bool need_finalized_vars)
{
  auto vars = obj->getMooseVariableDependencies();
  for (auto var : vars)
  {
    auto key = "variable_" + var->name();
    auto live_key = "live_variable_" + var->name();
    auto finalized_key = "finalized_variable_" + var->name();
    bool have_finalized = gd.named_objects.count(finalized_key) > 0 &&
                          gd.named_objects[finalized_key].count(n->loopType()) > 0;
    bool have_live =
        gd.named_objects.count(live_key) > 0 && gd.named_objects[live_key].count(n->loopType()) > 0;
    // depend on the live version of the var if there is one and it deoesn't depend on n
    if (need_finalized_vars && have_finalized &&
        !gd.named_objects[finalized_key][n->loopType()]->dependsOn(n))
      n->needs(gd.named_objects[finalized_key][n->loopType()]);
    else if (need_live_vars && have_live &&
             !gd.named_objects[live_key][n->loopType()]->dependsOn(n))
      n->needs(gd.named_objects[live_key][n->loopType()]);
    else
    {
      // fallback to depending on the non-live/lagged version of the var
      if (gd.named_objects.count(key) > 0 && gd.named_objects[key].count(n->loopType()) > 0)
      {
        auto var_node = gd.named_objects[key][n->loopType()];
        n->needs(var_node);
      }
      else
        mooseError("object's needed variable has no dag node");
    }
  }
}

void
addMatDeps(GraphData & gd, MaterialPropertyInterface * obj, dag::Node * n)
{
  auto & props = obj->getMatPropDependencies();
  for (auto prop : props)
  {
    mooseAssert(gd.named_mat_props.count(prop) > 0 &&
                    gd.named_mat_props[prop].count(n->loopType()) > 0,
                "object's needed material has no dag node");
    auto mat = gd.named_mat_props[prop][n->loopType()];
    n->needs(mat);
  }
}

// removes tail node from g (if it exists in g) if and only if head node is
// not present in g.  This helps with things like setup+teardown nodes where
// we eneded up not having any nodes that depended on setup so it is missing,
// but something still (unnecessarily) depended on the tail node.  We want to
// prune these nodes away.
void
pruneBrokenSandwiches(dag::Subgraph & g, dag::Node * head, dag::Node * tail)
{
  if (!g.contains(head))
  {
    g.remove(tail);
    tail->clearDeps();
  }
}

// convert the given (thread copies) of the given variable to a dag node of
// the given loop type.  We call this function several times for a given
// variable in order to create one dag node for every loop type this variable
// could be used in - although each of these copies is bound to the same moose
// object/run-func calls.  If live is true, then the dag node name will be the object name
// prefixed with "live_varible_" instead of just "variable_".  This means that
// the variable represents one that is "current" or "updated" for the current
// iteration - i.e. it has been solved on the current iteration prior to this
// variable node being calculated.  This is relevant for auxiliary variables
// that may have dependers who expect to see the non-lagged,
// computed-this-iteration value.  This could also be relevant for dampers and
// regular nonlinear variables if we ever get those implemented.
dag::Node *
convertVar(FEProblemBase &,
           GraphData & gd,
           std::vector<MooseVariableFieldBase *> & vars,
           dag::LoopType type,
           bool live = false, bool finalized = false)
{
  std::vector<dag::Node *> objs;
  auto nonlivename = "variable_" + vars[0]->name();
  auto objname = nonlivename;
  if (live)
  {
    objname = "live_" + nonlivename;
    if (finalized)
      objname = "finalized_" + nonlivename;
  }

  auto obj = gd.graph.create(objname, false, false, type);

  if (live && finalized)
  {
    // create an intermediate aux solve node that the variable node depends
    // on.  The auxsolve node also depends on the aux kernel for this variable.
    // This intermediate node has a None type and depends on a full aux
    // solution finalizing operation that is also None type - so it can run
    // together with the full aux system solve node in between the actual variable node and the
    // kernel. This 3 layer approach allows multiple loops depending on various aux variables to
    // share a single aux system solve operation that runs once. There may be multiple aux system
    // solve ops for the entire simulation per iteration, but they will be shared as much as
    // possible (hopefully) with this approach.
    auto solvname = "auxsolve_" + nonlivename;
    if (gd.named_objects.count(solvname) == 0)
    {
      auto solv = gd.graph.create(solvname, true, false, NoneType);
      gd.named_objects[solv->name()][solv->loopType()] = solv;
      solv->preserve();
    }
    auto solv = gd.named_objects[solvname][NoneType];
    obj->needs(solv);
    solv->needs(gd.aux_soln_finalize);
  }

  obj->setRunFunc([vars](const MeshLocation & loc, THREAD_ID tid) {
    auto var = vars[tid];
    var->clearDofIndices();
    var->prepare();
    if (loc.type.category == dag::LoopCategory::Elemental_onElem)
      var->computeElemValues();
    else if (loc.type.category == dag::LoopCategory::Elemental_onBoundary)
      var->computeElemValuesFace();
    else if (loc.type.category == dag::LoopCategory::Nodal)
    {
      var->reinitNode();
      if (var->isNodalDefined())
        var->computeNodalValues();
    }
    else if (loc.type.category == dag::LoopCategory::Nodal_onBoundary)
    {
      var->reinitNode();
      if (var->isNodalDefined())
        var->computeNodalValues();
    }
    else
      mooseError("unsupported loop type for variable node");
  });

  if (type.category == dag::LoopCategory::Elemental_onElem)
    obj->needs(gd.elem_setup[type.block]);
  else if (type.category == dag::LoopCategory::Elemental_onBoundary)
    obj->needs(gd.side_setup[type.block]);
  else if (type.category == dag::LoopCategory::Nodal_onBoundary)
    obj->needs(gd.side_node_setup[type.block]);
  else if (type.category == dag::LoopCategory::Nodal)
    obj->needs(gd.node_setup[type.block]);
  else
    mooseError("unsupported loop type for variable node");
  return obj;
}

dag::Node *
convertAuxKernel(FEProblemBase & fe,
                 GraphData & gd,
                 std::vector<AuxKernel *> & kernels,
                 dag::LoopType type)
{
  auto obj_name = "auxkernel_" + kernels[0]->name();
  auto obj = gd.graph.create(obj_name, false, false, type);
  obj->setRunFunc([&fe, kernels](const MeshLocation &, THREAD_ID tid) {
    kernels[tid]->compute();
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      kernels[tid]->variable().insert(fe.getAuxiliarySystem().solution());
    }
  });

  // make all loop type incarnations of the kernel's primary variable's "live"
  // dag node (on the current block) depend on this aux kernel. Do this before
  // we add variable dependencies - so that when addVarDeps tries to make this
  // aux kernel depend on its own variable (sloppy moose world) - it triggers
  // a cyclical dependency and can be skipped.
  //
  // It's also possible that the looptype blocks mismatch because e.g. a nodal
  // aux variable may depend on e.g. a boundary restricted aux kernel - if
  // this is the case, then we still need+want the var->aux dependency in
  // place - and in fact we want the variable to depend on every aux kernel
  // for that for the corresponding aux var - not just the block-matching one
  // (since there isn't really a block matching one).  So that is why we still
  // add deps below when the looptype categories mismatch (even if the blocks
  // also mismatch).  And in the case of live variables, we sort of want to
  // "promote" live variables to also depend on the auxsolve node like
  // finalized variables do in the case where the loop category types do
  // mismatch.
  for (auto & entry : gd.named_objects["live_variable_" + kernels[0]->variable().name()])
    if (entry.first.block == type.block && entry.first.category == type.category)
      entry.second->needs(obj);
    else if (entry.first.category != type.category)
    {
      entry.second->needs(obj);
      entry.second->needs(
          gd.named_objects["auxsolve_variable_" + kernels[0]->variable().name()][NoneType]);
    }
  for (auto & entry : gd.named_objects["finalized_variable_" + kernels[0]->variable().name()])
    if (entry.first.block == type.block || entry.first.category != type.category)
      entry.second->needs(obj);
  gd.named_objects["auxsolve_variable_" + kernels[0]->variable().name()][NoneType]->needs(obj);

  // TODO: we actually want live vars for variables of the same loop type and
  // finalized variables for variables of a different loop type.  Figure out
  // how to distinguish these cases.
  addVarDeps(gd, kernels[0], obj, true, false);
  addMatDeps(gd, kernels[0], obj);

  if (type.category == dag::LoopCategory::Nodal_onBoundary)
    obj->needs(gd.side_node_setup[type.block]);
  else if (type.category == dag::LoopCategory::Nodal)
  {
    obj->needs(gd.node_setup[type.block]);
    gd.node_teardown[type.block]->needs(obj);
  }
  else
    mooseError("unsupported loop type for auxvariable node");

  return obj;
}

// convert the given (thread copies) of the given material to a dag node of
// the given loop type.  We call this function several times for a given
// material in order to create one dag node for every loop type this material
// could be used in - although each of these copies is bound to the same moose
// object/run-func calls.
//
// We create a map of every material property the material defines to the dag node
// we create for this material.
dag::Node *
convertMat(FEProblemBase &, GraphData & gd, std::vector<Material *> & mats, dag::LoopType type)
{
  auto objname = "material_" + mats[0]->name();
  auto obj = gd.graph.create(objname, false, false, type);
  obj->setRunFunc([mats](const MeshLocation &, THREAD_ID tid) { mats[tid]->computeProperties(); });

  auto mat = mats[0];
  auto & props = mat->getSuppliedPropIDs();
  for (auto prop : props)
    gd.named_mat_props[prop][type] = obj;

  // TODO: what if a kernel depends on this material and needs the "live"
  // version of a variable for which this material depends on the "lagged"
  // version? - this will cause both the live and lagged dag nodes to be in
  // the loop's objects - and so the variable will be computed twice - can we
  // make those variable->compute... calls idempotent to resolve this?
  addVarDeps(gd, mats[0], obj, false, false);

  // we can't add material->material dependencies in here because they might not have
  // been converted/created yet.

  if (type.category == dag::LoopCategory::Elemental_onElem)
    obj->needs(gd.elem_setup[type.block]);
  else if (type.category == dag::LoopCategory::Elemental_onBoundary)
    obj->needs(gd.side_setup[type.block]);
  else if (type.category == dag::LoopCategory::Nodal_onBoundary)
    obj->needs(gd.side_node_setup[type.block]);
  else
    mooseError("unsupported loop type for variable node");

  return obj;
}

dag::Node *
convertKernel(FEProblemBase &,
              GraphData & gd,
              std::vector<KernelBase *> & kernels,
              SubdomainID block)
{
  auto objname = "kernel_" + kernels[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
  obj->setRunFunc(
      [kernels](const MeshLocation &, THREAD_ID tid) { kernels[tid]->computeResidual(); });

  addVarDeps(gd, kernels[0], obj, true, true);
  addMatDeps(gd, kernels[0], obj);

  if (dynamic_cast<TimeKernel *>(kernels[0]))
    obj->needs(gd.time_derivative);

  obj->needs(gd.elem_setup[block]);
  gd.elem_teardown[block]->needs(obj);
  obj->needs(gd.residual_setup);
  gd.residual_teardown->needs(obj);
  return obj;
}

dag::Node *
convertBC(FEProblemBase &,
          GraphData & gd,
          std::vector<IntegratedBCBase *> & bcs,
          BoundaryID boundary)
{
  auto objname = "bc_" + bcs[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
  obj->setRunFunc([bcs](const MeshLocation &, THREAD_ID tid) {
    if (bcs[tid]->shouldApply())
      bcs[tid]->computeResidual();
  });

  addVarDeps(gd, bcs[0], obj, true, true);
  addMatDeps(gd, bcs[0], obj);

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
  auto objname = "bc_" + bcs[0]->name();
  auto obj = gd.graph.create(
      objname, false, false, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
  obj->setRunFunc([bcs](const MeshLocation & /*loc*/, THREAD_ID tid) {
    if (bcs[tid]->shouldApply())
      bcs[tid]->computeResidual();
  });

  addVarDeps(gd, bcs[0], obj, true, true);

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
  // create variable nodes
  {
    auto field_vars = fe.getNonlinearSystemBase().getVariables(tmp_tid);
    for (auto v : field_vars)
    {
      std::vector<MooseVariableFieldBase *> var_copies;
      for (unsigned int i = 0; i < n_threads; i++)
        var_copies.push_back(&fe.getNonlinearSystemBase().getVariable(i, v->name()));

      dag::Node * obj = nullptr;
      // here we loop over all mesh subdomains instead of just the ones the
      // variable is defined on because then we don't have to worry about
      // handling the "special" subdomain ID values (ANY, etc.) - and
      // evaluation only occurs when we want it to if a particular block
      // restricted depender object grabs the var on this block.
      for (auto block : fe.mesh().meshSubdomains())
      {
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal, block));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
      }
      for (auto boundary : fe.mesh().meshBoundaryIds())
      {
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
      }
    }
  }

  // create aux variable nodes
  {
    auto field_vars = fe.getAuxiliarySystem().getVariables(tmp_tid);
    for (auto v : field_vars)
    {
      std::vector<MooseVariableFieldBase *> var_copies;
      for (unsigned int i = 0; i < n_threads; i++)
        var_copies.push_back(&fe.getAuxiliarySystem().getVariable(i, v->name()));

      dag::Node * obj = nullptr;
      // here we loop over all mesh subdomains instead of just the ones the
      // variable is defined on because then we don't have to worry about
      // handling the "special" subdomain ID values (ANY, etc.) - and
      // evaluation only occurs when we want it to if a particular block
      // restricted depender object grabs the var on this block.
      for (auto block : fe.mesh().meshSubdomains())
      {
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal, block));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        // we also create "live" versions of each aux var node - that aux
        // kernels will add dependencies from - on themselves.  This allows
        // objects to depend on either lagged or current values of aux variables.
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onElem, block), true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal, block), true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onElem, block), true, true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal, block), true, true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
      }
      for (auto boundary : fe.mesh().meshBoundaryIds())
      {
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        // we also create "live" versions of each aux var node - that aux
        // kernels will add dependencies from - on themselves.  This allows
        // objects to depend on either lagged or current values of aux variables.
        obj = convertVar(fe,
                         gd,
                         var_copies,
                         dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary),
                         true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary), true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(fe,
                         gd,
                         var_copies,
                         dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary),
                         true, true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
        obj = convertVar(
            fe, gd, var_copies, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary), true, true);
        gd.named_objects[obj->name()][obj->loopType()] = obj;
      }
    }
  }

  // create material nodes
  std::vector<std::pair<MaterialPropertyInterface *, dag::Node *>> material_node_pairs;
  auto & mat_warehouse = fe.getMaterialWarehouse();
  for (auto mat : mat_warehouse.getActiveObjects(tmp_tid))
  {
    auto m = dynamic_cast<Material *>(mat.get());
    std::vector<Material *> mat_copies;
    for (unsigned int i = 0; i < n_threads; i++)
      mat_copies.push_back(dynamic_cast<Material *>(mat_warehouse.getObject(m->name(), i).get()));
    for (auto block : fe.mesh().meshSubdomains())
    {
      auto obj =
          convertMat(fe, gd, mat_copies, dag::LoopType(dag::LoopCategory::Elemental_onElem, block));
      material_node_pairs.emplace_back(m, obj);
    }
    for (auto boundary : fe.mesh().meshBoundaryIds())
    {
      auto obj = convertMat(
          fe, gd, mat_copies, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
      material_node_pairs.emplace_back(m, obj);
    }
  }
  for (auto mat : mat_warehouse[Moose::FACE_MATERIAL_DATA].getActiveObjects(tmp_tid))
  {
    auto m = dynamic_cast<Material *>(mat.get());
    std::vector<Material *> mat_copies;
    for (unsigned int i = 0; i < n_threads; i++)
      mat_copies.push_back(dynamic_cast<Material *>(
          mat_warehouse[Moose::FACE_MATERIAL_DATA].getObject(m->name(), i).get()));
    for (auto boundary : fe.mesh().meshBoundaryIds())
    {
      auto obj = convertMat(
          fe, gd, mat_copies, dag::LoopType(dag::LoopCategory::Elemental_onBoundary, boundary));
      material_node_pairs.emplace_back(m, obj);
    }
  }
  // add material->material dependencies *after* we've created all material dag nodes
  for (auto & pair : material_node_pairs)
    addMatDeps(gd, pair.first, pair.second);

  // create kernels and bcs
  auto & auxkern_warehouse = fe.getAuxiliarySystem().getNodalKernelWarehouse();
  auto & kern_warehouse =
      fe.getNonlinearSystemBase().getKernelWarehouse().getVectorTagsObjectWarehouse(tags, tmp_tid);
  auto & ibc_warehouse =
      fe.getNonlinearSystemBase().getIntegratedBCWarehouse().getVectorTagsObjectWarehouse(tags,
                                                                                          tmp_tid);
  auto & nbc_warehouse =
      fe.getNonlinearSystemBase().getNodalBCWarehouse().getVectorTagsObjectWarehouse(tags, tmp_tid);

  for (auto block : fe.mesh().meshSubdomains())
  {
    if (kern_warehouse.hasActiveBlockObjects(block, tmp_tid))
    {
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
    if (auxkern_warehouse.hasActiveBlockObjects(block, tmp_tid))
    {
      auto & kernels = auxkern_warehouse.getBlockObjects(block, tmp_tid);
      for (auto kern : kernels)
      {
        auto k = kern.get();
        std::vector<AuxKernel *> kernel_copies;
        for (unsigned int i = 0; i < n_threads; i++)
          kernel_copies.push_back(auxkern_warehouse.getObject(k->name(), i).get());
        convertAuxKernel(fe, gd, kernel_copies, dag::LoopType(dag::LoopCategory::Nodal, block));
      }
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
    if (auxkern_warehouse.hasActiveBoundaryObjects(boundary, tmp_tid))
    {
      auto & kernels = auxkern_warehouse.getBoundaryObjects(boundary, tmp_tid);
      for (auto kern : kernels)
      {
        auto k = kern.get();
        std::vector<AuxKernel *> kernel_copies;
        for (unsigned int i = 0; i < n_threads; i++)
          kernel_copies.push_back(auxkern_warehouse.getObject(k->name(), i).get());
        convertAuxKernel(
            fe, gd, kernel_copies, dag::LoopType(dag::LoopCategory::Nodal_onBoundary, boundary));
      }
    }
  }

  // prune un-needed stuff
  dag::pruneUp(gd.graph);
  for (auto block : fe.mesh().meshSubdomains())
    pruneBrokenSandwiches(gd.graph, gd.elem_setup[block], gd.elem_teardown[block]);
  for (auto boundary : fe.mesh().meshBoundaryIds())
    pruneBrokenSandwiches(gd.graph, gd.side_setup[boundary], gd.side_teardown[boundary]);

  std::cout << dag::dotGraph(gd.graph);

  gd.partitions = dag::computePartitions(gd.graph, true);

  gd.objs = computeLoops(gd.partitions);
  for (auto nodes : gd.objs)
    gd.loop_type.push_back(nodes[0][0]->loopType());
  dag::printLoops(gd.objs);

  // show all the given subgraphs on a single graph.
  // std::cout << dag::dotGraphMerged(gd.partitions);
}

} // namespace translation
