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

#ifndef MOOSEVARIABLEFAST_H
#define MOOSEVARIABLEFAST_H

#include "MooseVariable.h"

#include "MooseVariable.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "MooseVariableFast.h"

// libMesh
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/dense_vector.h"

/**
 * Class for stuff related to variables
 *
 * Each variable can compute nodal or elemental (at QPs) values.
 */
template <
  bool is_transient,

  bool need_u_old,
  bool need_u_older,
  bool need_u_previous_nl,

  bool need_grad_old,
  bool need_grad_older,
  bool need_grad_previous_nl,

  bool need_second,
  bool need_second_old,
  bool need_second_older,
  bool need_second_previous_nl,

  bool need_u_old_neighbor,
  bool need_u_older_neighbor,
  bool need_u_previous_nl_neighbor,

  bool need_grad_old_neighbor,
  bool need_grad_older_neighbor,
  bool need_grad_previous_nl_neighbor,

  bool need_second_neighbor,
  bool need_second_old_neighbor,
  bool need_second_older_neighbor,
  bool need_second_previous_nl_neighbor,

  bool need_nodal_u,
  bool need_nodal_u_old,
  bool need_nodal_u_older,
  bool need_nodal_u_previous_nl,
  bool need_nodal_u_dot,

  bool need_nodal_u_neighbor,
  bool need_nodal_u_old_neighbor,
  bool need_nodal_u_older_neighbor,
  bool need_nodal_u_previous_nl_neighbor,
  bool need_nodal_u_dot_neighbor
>
class MooseVariableFast : public MooseVariable
{
public:
  MooseVariableFast(MooseVariable * ref, Assembly& assem) :
    MooseVariable::MooseVariable(ref->number(), ref->feType(), ref->sys(), assem, ref->kind())
  {
  };

  virtual void computeElemValues() override
  {
    unsigned int nqp = _qrule->n_points();

    _u.resize(nqp);
    _grad_u.resize(nqp);

    if (need_second)
      _second_u.resize(nqp);

    if (need_u_previous_nl)
      _u_previous_nl.resize(nqp);

    if (need_grad_previous_nl)
      _grad_u_previous_nl.resize(nqp);

    if (need_second_previous_nl)
      _second_u_previous_nl.resize(nqp);

    if (is_transient)
    {
      _u_dot.resize(nqp);
      _du_dot_du.resize(nqp);

      if (need_u_old)
        _u_old.resize(nqp);

      if (need_u_older)
        _u_older.resize(nqp);

      if (need_grad_old)
        _grad_u_old.resize(nqp);

      if (need_grad_older)
        _grad_u_older.resize(nqp);

      if (need_second_old)
        _second_u_old.resize(nqp);

      if (need_second_older)
        _second_u_older.resize(nqp);
    }

    for (unsigned int i = 0; i < nqp; ++i)
    {
      _u[i] = 0;
      _grad_u[i] = 0;

      if (need_second)
        _second_u[i] = 0;

      if (need_u_previous_nl)
        _u_previous_nl[i] = 0;

      if (need_grad_previous_nl)
        _grad_u_previous_nl[i] = 0;

      if (need_second_previous_nl)
        _second_u_previous_nl[i] = 0;

      if (is_transient)
      {
        _u_dot[i] = 0;
        _du_dot_du[i] = 0;

        if (need_u_old)
          _u_old[i] = 0;

        if (need_u_older)
          _u_older[i] = 0;

        if (need_grad_old)
          _grad_u_old[i] = 0;

        if (need_grad_older)
          _grad_u_older[i] = 0;

        if (need_second_old)
          _second_u_old[i] = 0;

        if (need_second_older)
          _second_u_older[i] = 0;
      }
    }

    unsigned int num_dofs = _dof_indices.size();

    if (need_nodal_u)
      _nodal_u.resize(num_dofs);

    if (need_nodal_u_previous_nl)
      _nodal_u_previous_nl.resize(num_dofs);

    if (is_transient)
    {
      if (need_nodal_u_old)
        _nodal_u_old.resize(num_dofs);
      if (need_nodal_u_older)
        _nodal_u_older.resize(num_dofs);
      if (need_nodal_u_dot)
        _nodal_u_dot.resize(num_dofs);
    }

    const NumericVector<Real> & current_solution = *_sys.currentSolution();
    const NumericVector<Real> & solution_old     = _sys.solutionOld();
    const NumericVector<Real> & solution_older   = _sys.solutionOlder();
    const NumericVector<Real> * solution_prev_nl = _sys.solutionPreviousNewton();
    const NumericVector<Real> & u_dot            = _sys.solutionUDot();
    const Real & du_dot_du                       = _sys.duDotDu();

    dof_id_type idx = 0;
    Real soln_local = 0;
    Real soln_old_local = 0;
    Real soln_older_local = 0;
    Real soln_previous_nl_local = 0;
    Real u_dot_local = 0;

    Real phi_local = 0;
    const RealGradient * dphi_qp = NULL;
    const RealTensor * d2phi_local = NULL;

    RealGradient * grad_u_qp = NULL;

    RealGradient * grad_u_old_qp = NULL;
    RealGradient * grad_u_older_qp = NULL;
    RealGradient * grad_u_previous_nl_qp = NULL;

    RealTensor * second_u_qp = NULL;

    RealTensor * second_u_old_qp = NULL;
    RealTensor * second_u_older_qp = NULL;
    RealTensor * second_u_previous_nl_qp = NULL;

    for (unsigned int i=0; i < num_dofs; i++)
    {
      idx = _dof_indices[i];
      soln_local = current_solution(idx);

      if (need_nodal_u)
        _nodal_u[i] = soln_local;

      if (need_u_previous_nl || need_grad_previous_nl || need_second_previous_nl || need_nodal_u_previous_nl)
        soln_previous_nl_local = (*solution_prev_nl)(idx);

      if (need_nodal_u_previous_nl)
        _nodal_u_previous_nl[i] = soln_previous_nl_local;

      if (is_transient)
      {
        if (need_u_old || need_grad_old || need_second_old || need_nodal_u_old)
          soln_old_local = solution_old(idx);

        if (need_u_older || need_grad_older || need_second_older || need_nodal_u_older)
          soln_older_local = solution_older(idx);

        if (need_nodal_u_old)
          _nodal_u_old[i] = soln_old_local;
        if (need_nodal_u_older)
          _nodal_u_older[i] = soln_older_local;

        u_dot_local        = u_dot(idx);
        if (need_nodal_u_dot)
          _nodal_u_dot[i] = u_dot_local;
      }

      for (unsigned int qp=0; qp < nqp; qp++)
      {
        phi_local = _phi[i][qp];
        dphi_qp = &_grad_phi[i][qp];

        grad_u_qp = &_grad_u[qp];

        if (need_grad_previous_nl)
          grad_u_previous_nl_qp = &_grad_u_previous_nl[qp];

        if (is_transient)
        {
          if (need_grad_old)
            grad_u_old_qp = &_grad_u_old[qp];

          if (need_grad_older)
            grad_u_older_qp = &_grad_u_older[qp];
        }

        if (need_second || need_second_old || need_second_older || need_second_previous_nl)
        {
          d2phi_local = &(*_second_phi)[i][qp];

          if (need_second)
          {
            second_u_qp = &_second_u[qp];
            second_u_qp->add_scaled(*d2phi_local, soln_local);
          }

          if (need_second_previous_nl)
          {
            second_u_previous_nl_qp = &_second_u_previous_nl[qp];
            second_u_previous_nl_qp->add_scaled(*d2phi_local, soln_previous_nl_local);
          }

          if (is_transient)
          {
            if (need_second_old)
              second_u_old_qp = &_second_u_old[qp];

            if (need_second_older)
              second_u_older_qp = &_second_u_older[qp];
          }
        }

        _u[qp] += phi_local * soln_local;

        grad_u_qp->add_scaled(*dphi_qp, soln_local);

        if (need_u_previous_nl)
          _u_previous_nl[qp] += phi_local * soln_previous_nl_local;
        if (need_grad_previous_nl)
          grad_u_previous_nl_qp->add_scaled(*dphi_qp, soln_previous_nl_local);

        if (is_transient)
        {
          _u_dot[qp]        += phi_local * u_dot_local;
          _du_dot_du[qp]    = du_dot_du;

          if (need_u_old)
            _u_old[qp]        += phi_local * soln_old_local;

          if (need_u_older)
            _u_older[qp]      += phi_local * soln_older_local;

          if (need_grad_old)
            grad_u_old_qp->add_scaled(*dphi_qp, soln_old_local);

          if (need_grad_older)
            grad_u_older_qp->add_scaled(*dphi_qp, soln_older_local);

          if (need_second_old)
            second_u_old_qp->add_scaled(*d2phi_local, soln_old_local);

          if (need_second_older)
            second_u_older_qp->add_scaled(*d2phi_local, soln_older_local);
        }
      }
    }
  }
};


#endif /* MOOSEVARIABLEFAST_H */
