//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseVariableData.h"
#include "MooseArray.h"
#include "MooseTypes.h"

#include "libmesh/tensor_tools.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_n_tensor.h"
#include "libmesh/fe_type.h"
#include "DualRealOps.h"

#include <functional>
#include <vector>

class FaceInfo;

template <typename>
class MooseVariableFV;

template <typename OutputType>
class MooseVariableDataFV
{
public:
  using OutputGradient = typename MooseVariableData<OutputType>::OutputGradient;
  using OutputSecond = typename MooseVariableData<OutputType>::OutputSecond;
  using OutputDivergence = typename MooseVariableData<OutputType>::OutputDivergence;

  using FieldVariableValue = typename MooseVariableData<OutputType>::FieldVariableValue;
  using FieldVariableGradient = typename MooseVariableData<OutputType>::FieldVariableGradient;
  using FieldVariableSecond = typename MooseVariableData<OutputType>::FieldVariableSecond;
  using FieldVariableCurl = typename MooseVariableData<OutputType>::FieldVariableCurl;
  using FieldVariableDivergence = typename MooseVariableData<OutputType>::FieldVariableDivergence;

  using OutputShape = typename MooseVariableData<OutputType>::OutputShape;
  using OutputShapeGradient = typename MooseVariableData<OutputType>::OutputShapeGradient;
  using OutputShapeSecond = typename MooseVariableData<OutputType>::OutputShapeSecond;
  using OutputShapeDivergence = typename MooseVariableData<OutputType>::OutputShapeDivergence;

  using OutputData = typename MooseVariableData<OutputType>::OutputData;
  using DoFValue = typename MooseVariableData<OutputType>::DoFValue;

  MooseVariableDataFV(const MooseVariableFV<OutputType> & var,
                      const SystemBase & sys,
                      THREAD_ID tid,
                      Moose::ElementType element_type,
                      const Elem * const & elem);

  /**
   * Returns whether this data structure needs automatic differentiation calculations
   */
  bool needsAD() const { return _need_ad; }

  /**
   * Set the geometry type before calculating variables values
   * @param gm_type The type type of geometry; either Volume or Face
   */
  void setGeometry(Moose::GeometryType gm_type);

  //////////////// Heavy lifting computational routines //////////////////////////////

  /**
   * compute the variable values
   */
  void computeValuesFace(const FaceInfo & fi);

  /**
   * compute the variable values
   */
  void computeValues(bool force = false);

  /**
   * compute AD things
   */
  void computeAD(const unsigned int num_dofs, const unsigned int nqp);

  ///////////////////////////// Shape functions /////////////////////////////////////

  /**
   * The current element
   */
  const Elem * const & currentElem() const { return _elem; }

  /**
   * prepare the initial condition
   */
  void prepareIC();

  //////////////////////////////////// Solution getters /////////////////////////////////////

  /**
   * Local solution getter
   * @param state The state of the simulation: current, old, older, previous nl
   */
  const FieldVariableValue & sln(Moose::SolutionState state) const;

  /**
   * Local solution gradient getter
   * @param state The state of the simulation: current, old, older, previous nl
   */
  const FieldVariableGradient & gradSln(Moose::SolutionState state) const;

  /**
   * Local time derivative of solution gradient getter
   */
  const FieldVariableGradient & gradSlnDot() const;

  /**
   * Local second time derivative of solution gradient getter
   */
  const FieldVariableGradient & gradSlnDotDot() const;

  /**
   * Local solution second spatial derivative getter
   * @param state The state of the simulation: current, old, older, previous nl
   */
  const FieldVariableSecond & secondSln(Moose::SolutionState state) const;

  /**
   * Local solution curl getter
   * @param state The state of the simulation: current, old, older
   */
  const FieldVariableCurl & curlSln(Moose::SolutionState state) const;

  template <ComputeStage compute_stage>
  const typename VariableValueType<OutputType, compute_stage>::type & adSln() const
  {
    _need_ad = _need_ad_u = true;
    return _ad_u;
  }

  template <ComputeStage compute_stage>
  const typename VariableGradientType<OutputType, compute_stage>::type & adGradSln() const
  {
    _need_ad = _need_ad_grad_u = true;
    return _ad_grad_u;
  }

  template <ComputeStage compute_stage>
  const typename VariableSecondType<OutputType, compute_stage>::type & adSecondSln() const
  {
    _need_ad = _need_ad_second_u = true;
    return _ad_second_u;
  }

  template <ComputeStage compute_stage>
  const typename VariableValueType<OutputType, compute_stage>::type & adUDot() const
  {
    _need_ad = _need_ad_u_dot = true;
    return _ad_u_dot;
  }

  const FieldVariableValue & uDot() const;

  const FieldVariableValue & uDotDot() const;

  const FieldVariableValue & uDotOld() const;

  const FieldVariableValue & uDotDotOld() const;

  const VariableValue & duDotDu() const
  {
    _need_du_dot_du = true;
    return _du_dot_du;
  }

  const VariableValue & duDotDotDu() const
  {
    _need_du_dotdot_du = true;
    return _du_dotdot_du;
  }

  /**
   * Set local DOF values and evaluate the values on quadrature points
   */
  void setDofValues(const DenseVector<OutputData> & values);

  void setDofValue(const OutputData & value, unsigned int index);

  OutputData
  getElementalValue(const Elem * elem, Moose::SolutionState state, unsigned int idx = 0) const;

  ///////////////////////////// dof indices ///////////////////////////////////////////////

  virtual void getDofIndices(const Elem * elem, std::vector<dof_id_type> & dof_indices) const;
  const std::vector<dof_id_type> & dofIndices() const { return _dof_indices; }
  unsigned int numberOfDofs() const { return _dof_indices.size(); }
  void clearDofIndices() { _dof_indices.clear(); }

  /**
   * Set the current local DOF values to the input vector
   */
  void insert(NumericVector<Number> & residual);
  /**
   * Add the current local DOF values to the input vector
   */
  void add(NumericVector<Number> & residual);
  /**
   * Add passed in local DOF values to a solution vector
   */
  void addSolution(NumericVector<Number> & sol, const DenseVector<Number> & v) const;

  /////////////////////////// DoF value getters /////////////////////////////////////

  const DoFValue & dofValues() const;
  const DoFValue & dofValuesOld() const;
  const DoFValue & dofValuesOlder() const;
  const DoFValue & dofValuesPreviousNL() const;
  const DoFValue & dofValuesDot() const;
  const DoFValue & dofValuesDotOld() const;
  const DoFValue & dofValuesDotDot() const;
  const DoFValue & dofValuesDotDotOld() const;
  const MooseArray<Number> & dofValuesDuDotDu() const;
  const MooseArray<Number> & dofValuesDuDotDotDu() const;

  /**
   * Return the AD dof values
   */
  template <ComputeStage compute_stage>
  const MooseArray<typename Moose::RealType<compute_stage>::type> & adDofValues() const;

  /////////////////////////////// Increment stuff ///////////////////////////////////////

  /**
   * Increment getter
   * @return The increment
   */
  const FieldVariableValue & increment() const { return _increment; }

  /**
   * Compute and store incremental change in solution at QPs based on increment_vec
   */
  void computeIncrementAtQps(const NumericVector<Number> & increment_vec);

  /////////////////////////////// Tags ///////////////////////////////////////////////////

  const FieldVariableValue & vectorTagValue(TagID tag)
  {
    _need_vector_tag_u[tag] = true;
    return _vector_tag_u[tag];
  }
  const FieldVariableValue & matrixTagValue(TagID tag)
  {
    _need_matrix_tag_u[tag] = true;
    return _matrix_tag_u[tag];
  }

private:
  void initializeSolnVars();

  /**
   * Helper methods for assigning nodal values from their corresponding solution values (dof
   * values as they're referred to here in this class). These methods are only truly meaningful
   * for nodal basis families
   */
  void fetchDoFValues();
  void fetchADDoFValues();
  void zeroSizeDofValues();

  /// A const reference to the owning MooseVariableFE object
  const MooseVariableFV<OutputType> & _var;

  const FEType & _fe_type;

  const unsigned int _var_num;

  const SystemBase & _sys;

  const SubProblem & _subproblem;

  THREAD_ID _tid;

  const Assembly & _assembly;

  const DofMap & _dof_map;

  /// The element type this object is storing data for. This is either Element, Neighbor, or Lower
  Moose::ElementType _element_type;

  /// Number of components of the associated variable
  unsigned int _count;

  /// The dof indices for the current element
  std::vector<dof_id_type> _dof_indices;

  mutable std::vector<bool> _need_vector_tag_dof_u;
  mutable std::vector<bool> _need_matrix_tag_dof_u;

  // Dof values of tagged vectors
  std::vector<DoFValue> _vector_tags_dof_u;
  // Dof values of the diagonal of tagged matrices
  std::vector<DoFValue> _matrix_tags_dof_u;

  std::vector<FieldVariableValue> _vector_tag_u;
  mutable std::vector<bool> _need_vector_tag_u;
  std::vector<FieldVariableValue> _matrix_tag_u;
  mutable std::vector<bool> _need_matrix_tag_u;

  /// Continuity type of the variable
  FEContinuity _continuity;

  /// Increment in the variable used in dampers
  FieldVariableValue _increment;

  /// A zero AD variable
  const DualReal _ad_zero;

  /// u flags
  mutable bool _need_u_old;
  mutable bool _need_u_older;
  mutable bool _need_u_previous_nl;

  /// u dot flags
  mutable bool _need_u_dot;
  mutable bool _need_ad_u_dot;
  mutable bool _need_u_dotdot;
  mutable bool _need_u_dot_old;
  mutable bool _need_u_dotdot_old;
  mutable bool _need_du_dot_du;
  mutable bool _need_du_dotdot_du;

  /// gradient flags
  mutable bool _need_grad_old;
  mutable bool _need_grad_older;
  mutable bool _need_grad_previous_nl;

  /// gradient dot flags
  mutable bool _need_grad_dot;
  mutable bool _need_grad_dotdot;

  /// SolutionState second_u flags
  mutable bool _need_second;
  mutable bool _need_second_old;
  mutable bool _need_second_older;
  mutable bool _need_second_previous_nl;

  /// curl flags
  mutable bool _need_curl;
  mutable bool _need_curl_old;
  mutable bool _need_curl_older;

  /// AD flags
  mutable bool _need_ad;
  mutable bool _need_ad_u;
  mutable bool _need_ad_grad_u;
  mutable bool _need_ad_second_u;

  /// local solution flags
  mutable bool _need_dof_values;
  mutable bool _need_dof_values_old;
  mutable bool _need_dof_values_older;
  mutable bool _need_dof_values_previous_nl;
  mutable bool _need_dof_values_dot;
  mutable bool _need_dof_values_dotdot;
  mutable bool _need_dof_values_dot_old;
  mutable bool _need_dof_values_dotdot_old;
  mutable bool _need_dof_du_dot_du;
  mutable bool _need_dof_du_dotdot_du;

  /// local solution values
  DoFValue _dof_values;
  DoFValue _dof_values_old;
  DoFValue _dof_values_older;
  DoFValue _dof_values_previous_nl;

  /// nodal values of u_dot
  DoFValue _dof_values_dot;
  /// nodal values of u_dotdot
  DoFValue _dof_values_dotdot;
  /// nodal values of u_dot_old
  DoFValue _dof_values_dot_old;
  /// nodal values of u_dotdot_old
  DoFValue _dof_values_dotdot_old;
  /// nodal values of derivative of u_dot wrt u
  MooseArray<Number> _dof_du_dot_du;
  /// nodal values of derivative of u_dotdot wrt u
  MooseArray<Number> _dof_du_dotdot_du;

  /// u
  FieldVariableValue _u;
  FieldVariableValue _u_old;
  FieldVariableValue _u_older;
  FieldVariableValue _u_previous_nl;

  /// grad_u
  FieldVariableGradient _grad_u;
  FieldVariableGradient _grad_u_old;
  FieldVariableGradient _grad_u_older;
  FieldVariableGradient _grad_u_previous_nl;

  /// grad_u dots
  FieldVariableGradient _grad_u_dot;
  FieldVariableGradient _grad_u_dotdot;

  /// second_u
  FieldVariableSecond _second_u;
  FieldVariableSecond _second_u_old;
  FieldVariableSecond _second_u_older;
  FieldVariableSecond _second_u_previous_nl;

  /// curl_u
  FieldVariableCurl _curl_u;
  FieldVariableCurl _curl_u_old;
  FieldVariableCurl _curl_u_older;

  /// AD u
  typename VariableValueType<OutputShape, JACOBIAN>::type _ad_u;
  typename VariableGradientType<OutputShape, JACOBIAN>::type _ad_grad_u;
  typename VariableSecondType<OutputShape, JACOBIAN>::type _ad_second_u;
  MooseArray<DualReal> _ad_dof_values;
  MooseArray<DualReal> _ad_dofs_dot;
  typename VariableValueType<OutputShape, JACOBIAN>::type _ad_u_dot;

  // time derivatives

  /// u_dot (time derivative)
  FieldVariableValue _u_dot;

  /// u_dotdot (second time derivative)
  FieldVariableValue _u_dotdot, _u_dotdot_bak;

  /// u_dot_old (time derivative)
  FieldVariableValue _u_dot_old, _u_dot_old_bak;

  /// u_dotdot_old (second time derivative)
  FieldVariableValue _u_dotdot_old, _u_dotdot_old_bak;

  /// derivative of u_dot wrt u
  VariableValue _du_dot_du;

  /// derivative of u_dotdot wrt u
  VariableValue _du_dotdot_du, _du_dotdot_du_bak;

  /// Pointer to time integrator
  const TimeIntegrator * _time_integrator;

  /// The current elem. This has to be a reference because the current elem will be constantly
  /// changing. If we initialized this to point to one elem, then in the next calculation we would
  /// be pointing to the wrong place!
  const Elem * const & _elem;

  /// Whether this variable is being calculated on a displaced system
  const bool _displaced;
};

/////////////////////// General template definitions //////////////////////////////////////

template <typename OutputType>
template <ComputeStage compute_stage>
const MooseArray<typename Moose::RealType<compute_stage>::type> &
MooseVariableDataFV<OutputType>::adDofValues() const
{
  _need_ad = true;
  return _ad_dof_values;
}

////////////////////////// Forward declaration of fully specialized templates //////////////////

template <>
template <>
const VariableValue & MooseVariableDataFV<Real>::adSln<RESIDUAL>() const;

template <>
template <>
const VariableGradient & MooseVariableDataFV<Real>::adGradSln<RESIDUAL>() const;

template <>
template <>
const VariableSecond & MooseVariableDataFV<Real>::adSecondSln<RESIDUAL>() const;

template <>
template <>
const VariableValue & MooseVariableDataFV<Real>::adUDot<RESIDUAL>() const;

template <>
template <>
const VectorVariableValue & MooseVariableDataFV<RealVectorValue>::adSln<RESIDUAL>() const;

template <>
template <>
const VectorVariableGradient & MooseVariableDataFV<RealVectorValue>::adGradSln<RESIDUAL>() const;

template <>
template <>
const VectorVariableSecond & MooseVariableDataFV<RealVectorValue>::adSecondSln<RESIDUAL>() const;

template <>
template <>
const VectorVariableValue & MooseVariableDataFV<RealVectorValue>::adUDot<RESIDUAL>() const;

template <>
template <>
const MooseArray<Real> & MooseVariableDataFV<Real>::adDofValues<RESIDUAL>() const;

template <>
template <>
const MooseArray<Real> & MooseVariableDataFV<RealVectorValue>::adDofValues<RESIDUAL>() const;

