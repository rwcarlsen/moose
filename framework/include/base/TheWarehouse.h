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
  System,
  Enabled,
  Tag,       // multiple
  Boundary,  // multiple
  Subdomain, // multiple
  ExecOn,    // multiple
};

struct Attribute
{
  AttributeId id;
  int value;
  std::string strvalue;
  inline bool operator==(const Attribute & other) const
  {
    return id == other.id && value == other.value && strvalue == other.strvalue;
  }
};

class Warehouse
{
public:
  Warehouse();

  void add(std::unique_ptr<MooseObject> obj);

  // prepares a query and returns an associated query_id (i.e. for use with the query function.
  int prepare(const std::vector<Attribute> & conds);

  const std::vector<MooseObject *> & query(int query_id);

private:
  void readAttribs(MooseObject * std::vector<Attribute> & attribs)

      std::unique_ptr<Storage> _store;
  std::vector<std::unique_ptr<MooseObject>> _objects;

  std::vector<std::vector<MooseObject *>> _obj_cache;
  std::vector<std::vector<Attribute>> _query_cache;
  std::vector<bool> _query_dirty;
};

#endif // THEWAREHOUSE_H
