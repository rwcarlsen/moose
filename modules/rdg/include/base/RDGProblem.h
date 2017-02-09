/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef RDGPROBLEM_H
#define RDGPROBLEM_H

#include "FEProblemBase.h"

class RDGProblem;
class RDGSystem;

template<>
InputParameters validParams<RDGProblem>();


class RDGProblem : public FEProblemBase
{
public:
  RDGProblem(const InputParameters & parameters);
  virtual ~RDGProblem();

  virtual void solve() override;

  virtual void addUserObject(std::string user_object_name, const std::string & name, InputParameters parameters) override;

  RDGSystem * _rdg_sys;
};

#endif
