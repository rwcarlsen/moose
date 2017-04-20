/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef MOOSEVARIABLE_H
#define MOOSEVARIABLE_H

#include "MooseTypes.h"
#include "MooseVariableBase.h"

// Forward declarations
class Assembly;
class SubProblem;
class SystemBase;

// libMesh forward declarations
namespace libMesh
{
class Elem;
class Node;
class QBase;
template <typename T>
class NumericVector;
template <typename T>
class DenseVector;
}

/**
 * Class for stuff related to variables
 *
 * Each variable can compute nodal or elemental (at QPs) values.
 */
class MooseVariable : public MooseVariableBase
{
public:
  MooseVariable(unsigned int var_num,
                const FEType & fe_type,
                SystemBase & sys,
                Assembly & assembly,
                Moose::VarKindType var_kind);
  virtual ~MooseVariable();

  /**
   * Clear out the dof indices.  We do this in case this variable is not going to be prepared at
   * all...
   */
  void clearDofIndices();

  void prepare();

  void prepareNeighbor();
  void prepareAux();
  void prepareIC();

  void reinitNode();
  void reinitNodeNeighbor();
  void reinitAux();
  void reinitAuxNeighbor();

  void reinitNodes(const std::vector<dof_id_type> & nodes);
  void reinitNodesNeighbor(const std::vector<dof_id_type> & nodes);

  const std::set<SubdomainID> & activeSubdomains();

  /**
   * Is the variable active on the subdomain?
   * @param subdomain The subdomain id in question
   * @return true if active on subdomain, false otherwise
   */
  bool activeOnSubdomain(SubdomainID subdomain) const;

  /**
   * Is this variable nodal
   * @return true if it nodal, otherwise false
   */
  virtual bool isNodal() const override;

  /**
   * Current element this variable is evaluated at
   */
  const Elem *& currentElem() { return _elem; }

  /**
   * Current side this variable is being evaluated on
   */
  unsigned int & currentSide() { return _current_side; }

  /**
   * Current neighboring element
   */
  const Elem *& neighbor() { return _neighbor; }

  /**
   * Whether or not this variable is computing any second derivatives.
   */
  bool computingSecond()
  {
    return haveNeed(Need::second) || haveNeed(Need::second_old) || haveNeed(Need::second_older) || haveNeed(Need::second_previous_nl);
  }

  const VariablePhiValue & phi();
  const VariablePhiGradient & gradPhi();
  const VariablePhiSecond & secondPhi();

  const VariablePhiValue & phiFace();
  const VariablePhiGradient & gradPhiFace();
  const VariablePhiSecond & secondPhiFace();

  const VariablePhiValue & phiNeighbor();
  const VariablePhiGradient & gradPhiNeighbor();
  const VariablePhiSecond & secondPhiNeighbor();

  const VariablePhiValue & phiFaceNeighbor();
  const VariablePhiGradient & gradPhiFaceNeighbor();
  const VariablePhiSecond & secondPhiFaceNeighbor();

  const MooseArray<Point> & normals() { return _normals; }

  // damping
  VariableValue & increment() { return _increment; }

  const VariableValue & sln() { return _u; }
  const VariableValue & slnOld()
  {
    need(Need::u_old);
    return _u_old;
  }
  const VariableValue & slnOlder()
  {
    need(Need::u_older);
    return _u_older;
  }
  const VariableValue & slnPreviousNL()
  {
    need(Need::u_previous_nl);
    return _u_previous_nl;
  }
  const VariableGradient & gradSln() { return _grad_u; }
  const VariableGradient & gradSlnOld()
  {
    need(Need::grad_old);
    return _grad_u_old;
  }
  const VariableGradient & gradSlnOlder()
  {
    need(Need::grad_older);
    return _grad_u_older;
  }
  const VariableGradient & gradSlnPreviousNL()
  {
    need(Need::grad_previous_nl);
    return _grad_u_previous_nl;
  }
  const VariableSecond & secondSln()
  {
    need(Need::second);
    secondPhi();
    secondPhiFace();
    return _second_u;
  }
  const VariableSecond & secondSlnOld()
  {
    need(Need::second_old);
    secondPhi();
    secondPhiFace();
    return _second_u_old;
  }
  const VariableSecond & secondSlnOlder()
  {
    need(Need::second_older);
    secondPhi();
    secondPhiFace();
    return _second_u_older;
  }
  const VariableSecond & secondSlnPreviousNL()
  {
    need(Need::second_previous_nl);
    secondPhi();
    secondPhiFace();
    return _second_u_previous_nl;
  }

  const VariableValue & uDot() { return _u_dot; }
  const VariableValue & duDotDu() { return _du_dot_du; }

  const Node *& node() { return _node; }
  dof_id_type & nodalDofIndex() { return _nodal_dof_index; }
  bool isNodalDefined() { return _is_defined; }
  const VariableValue & nodalSln() { return _nodal_u; }
  const VariableValue & nodalSlnOld() { return _nodal_u_old; }
  const VariableValue & nodalSlnOlder() { return _nodal_u_older; }
  const VariableValue & nodalSlnPreviousNL()
  {
    need(Need::nodal_u_previous_nl);
    return _nodal_u_previous_nl;
  }
  const VariableValue & nodalSlnDot() { return _nodal_u_dot; }
  const VariableValue & nodalSlnDuDotDu() { return _nodal_du_dot_du; }

  const VariableValue & nodalValue();
  const VariableValue & nodalValueOld();
  const VariableValue & nodalValueOlder();
  const VariableValue & nodalValuePreviousNL();
  const VariableValue & nodalValueDot();

  const VariableValue & nodalValueNeighbor();
  const VariableValue & nodalValueOldNeighbor();
  const VariableValue & nodalValueOlderNeighbor();
  const VariableValue & nodalValuePreviousNLNeighbor();
  const VariableValue & nodalValueDotNeighbor();

  const VariableValue & slnNeighbor() { return _u_neighbor; }
  const VariableValue & slnOldNeighbor()
  {
    need(Need::u_old_neighbor);
    return _u_old_neighbor;
  }
  const VariableValue & slnOlderNeighbor()
  {
    need(Need::u_older_neighbor);
    return _u_older_neighbor;
  }
  const VariableValue & slnPreviousNLNeighbor()
  {
    need(Need::u_previous_nl_neighbor);
    return _u_previous_nl_neighbor;
  }
  const VariableGradient & gradSlnNeighbor() { return _grad_u_neighbor; }
  const VariableGradient & gradSlnOldNeighbor()
  {
    need(Need::grad_old_neighbor);
    return _grad_u_old_neighbor;
  }
  const VariableGradient & gradSlnOlderNeighbor()
  {
    need(Need::grad_older_neighbor);
    return _grad_u_older_neighbor;
  }
  const VariableGradient & gradSlnPreviousNLNeighbor()
  {
    need(Need::grad_previous_nl_neighbor);
    return _grad_u_previous_nl_neighbor;
  }
  const VariableSecond & secondSlnNeighbor()
  {
    need(Need::second_neighbor);
    secondPhiFaceNeighbor();
    return _second_u_neighbor;
  }
  const VariableSecond & secondSlnOldNeighbor()
  {
    need(Need::second_old_neighbor);
    secondPhiFaceNeighbor();
    return _second_u_old_neighbor;
  }
  const VariableSecond & secondSlnOlderNeighbor()
  {
    need(Need::second_older_neighbor);
    secondPhiFaceNeighbor();
    return _second_u_older_neighbor;
  }
  const VariableSecond & secondSlnPreviousNLNeighbor()
  {
    need(Need::second_previous_nl_neighbor);
    secondPhiFaceNeighbor();
    return _second_u_previous_nl_neighbor;
  }

  const VariableValue & uDotNeighbor() { return _u_dot_neighbor; }
  const VariableValue & duDotDuNeighbor() { return _du_dot_du_neighbor; }

  const Node *& nodeNeighbor() { return _node_neighbor; }
  dof_id_type & nodalDofIndexNeighbor() { return _nodal_dof_index_neighbor; }
  bool isNodalNeighborDefined() { return _is_defined_neighbor; }
  const VariableValue & nodalSlnNeighbor() { return _nodal_u_neighbor; }
  const VariableValue & nodalSlnOldNeighbor() { return _nodal_u_old_neighbor; }
  const VariableValue & nodalSlnOlderNeighbor() { return _nodal_u_older_neighbor; }
  const VariableValue & nodalSlnPreviousNLNeighbor()
  {
    need(Need::nodal_u_previous_nl_neighbor);
    return _nodal_u_previous_nl_neighbor;
  }
  const VariableValue & nodalSlnDotNeighbor() { return _nodal_u_dot_neighbor; }
  const VariableValue & nodalSlnDuDotDuNeighbor() { return _nodal_du_dot_du_neighbor; }

  /**
   * Compute values at interior quadrature points
   * when this variable's elem dof i is perturbed
   * by h in the current solution. Perturbation h
   * is computed based on the current variable
   * value at i, its size controlled by scale.
   */
  void computePerturbedElemValues(unsigned i, Real scale, Real & h);

  /**
   * Restore the values the variable had before a call to
   * computePerturbedElemValues().
   */
  void restoreUnperturbedElemValues();

  /**
   * Compute values at interior quadrature points
   */
  void computeElemValues();
  /**
   * Compute values at facial quadrature points
   */
  void computeElemValuesFace();
  /**
   * Compute values at facial quadrature points for the neighbor
   */
  void computeNeighborValuesFace();
  /**
   * Compute values at quadrature points for the neighbor
   */
  void computeNeighborValues();
  /**
   * Compute nodal values of this variable
   */
  void computeNodalValues();
  /**
   * Compute nodal values of this variable in the neighbor
   */
  void computeNodalNeighborValues();
  /**
   * Set the nodal value for this variable to keep everything up to date
   */
  void setNodalValue(Number value, unsigned int idx = 0);
  /**
   * Set values for this variable to keep everything up to date
   */
  void setNodalValue(const DenseVector<Number> & value);

  /**
   * Set the neighbor nodal value for this variable
   */
  void setNodalValueNeighbor(Number value);
  /**
   * Set the neighbor values for this variable
   */
  void setNodalValueNeighbor(const DenseVector<Number> & value);

  /**
   * Compute and store incremental change in solution at QPs based on increment_vec
   */
  void computeIncrementAtQps(const NumericVector<Number> & increment_vec);

  /**
   * Compute and store incremental change at the current node based on increment_vec
   */
  void computeIncrementAtNode(const NumericVector<Number> & increment_vec);

  /**
   * Get DOF indices for currently selected element
   * @return
   */
  std::vector<dof_id_type> & dofIndicesNeighbor() { return _dof_indices_neighbor; }

  unsigned int numberOfDofsNeighbor() { return _dof_indices_neighbor.size(); }

  void insert(NumericVector<Number> & residual);
  void add(NumericVector<Number> & residual);

  /**
   * Get the value of this variable at given node
   */
  Number getNodalValue(const Node & node);
  /**
   * Get the old value of this variable at given node
   */
  Number getNodalValueOld(const Node & node);
  /**
   * Get the t-2 value of this variable at given node
   */
  Number getNodalValueOlder(const Node & node);

  /**
   * Compute the variable value at a point on an element
   * @param elem The element we are computing on
   * @param phi Evaluated shape functions at a point
   * @return The variable value
   */
  Real getValue(const Elem * elem, const std::vector<std::vector<Real>> & phi) const;
  RealGradient getGradient(const Elem * elem,
                           const std::vector<std::vector<RealGradient>> & phi) const;

  /**
   * Retrieve the Elemental DOF
   * @param elem The element we are computing on
   * @return The variable value
   */
  Number getElementalValue(const Elem * elem, unsigned int idx = 0) const;

  /**
   * Whether or not this variable is actually using the shape function value.
   *
   * Currently hardcoded to true because we always compute the value.
   */
  bool usesPhi() { return true; }

  /**
   * Whether or not this variable is actually using the shape function gradient.
   *
   * Currently hardcoded to true because we always compute the value.
   */
  bool usesGradPhi() { return true; }

  /**
   * Whether or not this variable is actually using the shape function second derivative.
   */
  bool usesSecondPhi() { return haveNeed(Need::second) || haveNeed(Need::second_old) || haveNeed(Need::second_older); }

protected:
  /**
   * Get dof indices for the variable
   * @param elem Element whose DOFs we are requesting (input)
   * @param dof_indices DOF indices for the given element (output)
   */
  void getDofIndices(const Elem * elem, std::vector<dof_id_type> & dof_indices);

protected:
  /// Thread ID
  THREAD_ID _tid;

  /// Quadrature rule for interior
  QBase *& _qrule;
  /// Quadrature rule for the face
  QBase *& _qrule_face;
  /// Quadrature rule for the neighbor
  QBase *& _qrule_neighbor;

  /// current element
  const Elem *& _elem;
  /// the side of the current element (valid when doing face assembly)
  unsigned int & _current_side;

  /// neighboring element
  const Elem *& _neighbor;

  /// DOF indices (neighbor)
  std::vector<dof_id_type> _dof_indices_neighbor;

  // Shape function values, gradients. second derivatives
  const VariablePhiValue & _phi;
  const VariablePhiGradient & _grad_phi;
  const VariablePhiSecond * _second_phi;

  // Values, gradients and second derivatives of shape function on faces
  const VariablePhiValue & _phi_face;
  const VariablePhiGradient & _grad_phi_face;
  const VariablePhiSecond * _second_phi_face;

  // Values, gradients and second derivatives of shape function
  const VariablePhiValue & _phi_neighbor;
  const VariablePhiGradient & _grad_phi_neighbor;
  const VariablePhiSecond * _second_phi_neighbor;

  // Values, gradients and second derivatives of shape function on faces
  const VariablePhiValue & _phi_face_neighbor;
  const VariablePhiGradient & _grad_phi_face_neighbor;
  const VariablePhiSecond * _second_phi_face_neighbor;

  /// Normals at QPs on faces
  const MooseArray<Point> & _normals;

  VariableValue _u, _u_bak;
  VariableValue _u_old, _u_old_bak;
  VariableValue _u_older, _u_older_bak;
  VariableValue _u_previous_nl;
  VariableGradient _grad_u, _grad_u_bak;
  VariableGradient _grad_u_old, _grad_u_old_bak;
  VariableGradient _grad_u_older, _grad_u_older_bak;
  VariableGradient _grad_u_previous_nl;
  VariableSecond _second_u, _second_u_bak;
  VariableSecond _second_u_old, _second_u_old_bak;
  VariableSecond _second_u_older, _second_u_older_bak;
  VariableSecond _second_u_previous_nl;

  VariableValue _u_neighbor;
  VariableValue _u_old_neighbor;
  VariableValue _u_older_neighbor;
  VariableValue _u_previous_nl_neighbor;
  VariableGradient _grad_u_neighbor;
  VariableGradient _grad_u_old_neighbor;
  VariableGradient _grad_u_older_neighbor;
  VariableGradient _grad_u_previous_nl_neighbor;
  VariableSecond _second_u_neighbor;
  VariableSecond _second_u_old_neighbor;
  VariableSecond _second_u_older_neighbor;
  VariableSecond _second_u_previous_nl_neighbor;

  // time derivatives

  /// u_dot (time derivative)
  VariableValue _u_dot, _u_dot_bak;
  VariableValue _u_dot_neighbor, _u_dot_bak_neighbor;

  /// derivative of u_dot wrt u
  VariableValue _du_dot_du, _du_dot_du_bak;
  VariableValue _du_dot_du_neighbor, _du_dot_du_bak_neighbor;

  // nodal stuff

  /// If the variable is defined at the node (used in compute nodal values)
  bool _is_defined;
  /// If true, the nodal value gets inserted on calling insert()
  bool _has_nodal_value;
  bool _has_nodal_value_neighbor;
  const Node *& _node;
  dof_id_type _nodal_dof_index;
  VariableValue _nodal_u;
  VariableValue _nodal_u_old;
  VariableValue _nodal_u_older;
  VariableValue _nodal_u_previous_nl;

  /// nodal values of u_dot
  VariableValue _nodal_u_dot;
  /// nodal values of derivative of u_dot wrt u
  VariableValue _nodal_du_dot_du;

  /// If the variable is defined at the neighbor node (used in compute nodal values)
  bool _is_defined_neighbor;
  const Node *& _node_neighbor;
  dof_id_type _nodal_dof_index_neighbor;
  VariableValue _nodal_u_neighbor;
  VariableValue _nodal_u_old_neighbor;
  VariableValue _nodal_u_older_neighbor;
  VariableValue _nodal_u_previous_nl_neighbor;
  VariableValue _nodal_u_dot_neighbor;
  VariableValue _nodal_du_dot_du_neighbor;

  /// if variable is nodal
  bool _is_nodal;

  // damping
  VariableValue _increment;

  friend class NodeFaceConstraint;
  friend class ValueThresholdMarker;
  friend class ValueRangeMarker;

private:
  enum Need {
    none = 0,
    is_transient = 1 << 0,
    u_old = 1 << 1,
    u_older = 1 << 2,
    u_previous_nl = 1 << 3,
    grad_old = 1 << 4,
    grad_older = 1 << 5,
    grad_previous_nl = 1 << 6,
    second = 1 << 7,
    second_old = 1 << 8,
    second_older = 1 << 9,
    second_previous_nl = 1 << 10,
    u_old_neighbor = 1 << 11,
    u_older_neighbor = 1 << 12,
    u_previous_nl_neighbor = 1 << 13,
    grad_old_neighbor = 1 << 14,
    grad_older_neighbor = 1 << 15,
    grad_previous_nl_neighbor = 1 << 16,
    second_neighbor = 1 << 17,
    second_old_neighbor = 1 << 18,
    second_older_neighbor = 1 << 19,
    second_previous_nl_neighbor = 1 << 20,
    nodal_u = 1 << 21,
    nodal_u_old = 1 << 22,
    nodal_u_older = 1 << 23,
    nodal_u_previous_nl = 1 << 24,
    nodal_u_dot = 1 << 25,
    nodal_u_neighbor = 1 << 26,
    nodal_u_old_neighbor = 1 << 27,
    nodal_u_older_neighbor = 1 << 28,
    nodal_u_previous_nl_neighbor = 1 << 29,
    nodal_u_dot_neighbor = 1 << 30,
  };

  void resizeAll(unsigned int nqp, bool is_transient, unsigned int num_dofs);
  void resizeAllNeighbor(unsigned int nqp, bool is_transient, unsigned int num_dofs);
  bool tryFast();
  bool tryFastFace();

  inline bool haveNeed(Need n) const { return _needs & n; }
  void need(Need n) { _needs = static_cast<Need>(_needs | n); }

  template <int needs>
  void computeElemValuesFast();

  template <int needs>
  void computeElemValuesFaceFast();

  Need _needs = Need::none;
};

#endif /* MOOSEVARIABLE_H */
