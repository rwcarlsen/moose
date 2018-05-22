//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef THEWAREHOUSE_H
#define THEWAREHOUSE_H

#include <map>
#include <string>
#include <vector>

class MooseObject;
class Storage;

// attributes include:
//
//     * tag (multiple) few - 3ish
//     * system - order 50
//     * execute_on (multiple) 10 max
//     * thread_id - order 10
//     * boundary_id (multiple) 1000 per mesh, 1000 per object (use "all/any" optimization)
//     * subdomain_id (multiple) 10000 per mesh, 1000 per object (use "all/any" optimization)
//     * enabled
enum class AttributeId
{
  None,
  Thread,
  Type,
  Name,
  System,
  Variable,
  Enabled,
  Interfaces, // bitmask
  VectorTag,  // multiple
  MatrixTag,  // multiple
  Boundary,   // multiple
  Subdomain,  // multiple
  ExecOn,     // multiple
};

enum class Interfaces
{
  ElementUserObject,
  SideUserObject,
  InternalSideUserObject,
  NodalUserObject,
  GeneralUserObject,
  NonlocalKernel,
  NonlocalIntegratedBC,
  InternalSideIndicator,
  TransientMultiApp,
  MultiAppTransfer,
  Max, // This must be last
};

struct Attribute
{
  AttributeId id;
  int64_t value;
  std::string strvalue;
  inline bool operator==(const Attribute & other) const
  {
    return id == other.id && value == other.value && strvalue == other.strvalue;
  }
};

class TheWarehouse
{
public:
  TheWarehouse();
  ~TheWarehouse();

  void add(std::shared_ptr<MooseObject> obj, const std::string & system);
  void update(const MooseObject * obj, const std::string & system);

  // prepares a query and returns an associated query_id (i.e. for use with the query function.
  int prepare(const std::vector<Attribute> & conds);

  const std::vector<MooseObject *> & query(int query_id);

  template <typename T>
  void queryInto(int query_id, std::vector<T> & results)
  {
    auto & objs = query(query_id);
    results.resize(objs.size());
    for (unsigned int i = 0; i < objs.size(); i++)
    {
      auto obj = objs[i];
      assert(dynamic_cast<T *>(obj));
      results[i] = static_cast<T *>(obj);
    }
  }

private:
  void readAttribs(const MooseObject * obj,
                   const std::string & system,
                   std::vector<Attribute> & attribs);

  std::unique_ptr<Storage> _store;
  std::vector<std::shared_ptr<MooseObject>> _objects;

  std::vector<std::vector<MooseObject *>> _obj_cache;
  std::vector<std::vector<Attribute>> _query_cache;
  std::vector<bool> _query_dirty;
};

#endif // THEWAREHOUSE_H
