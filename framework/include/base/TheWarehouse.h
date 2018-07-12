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
#include <type_traits>
#include <unordered_map>
#include <iostream>

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
enum class AttributeId
{
  None,
  Thread,
  Name,
  System,
  Variable,
  Interfaces, // bitmask
  VectorTag,  // multiple
  MatrixTag,  // multiple
  Boundary,   // multiple
  Subdomain,  // multiple
  ExecOn,     // multiple
  // TODO: delete these two later - they are temporary hacks for dealing with inter-system
  // dependencies
  PreIC,
  PreAux,
};

enum class Interfaces
{
  ElementUserObject = 1 << 1,
  SideUserObject = 1 << 2,
  InternalSideUserObject = 1 << 3,
  NodalUserObject = 1 << 4,
  GeneralUserObject = 1 << 5,
  ShapeElementUserObject = 1 << 6,
  ShapeSideUserObject = 1 << 7,
  UserObject = 1 << 8,
  Postprocessor = 1 << 9,
  VectorPostprocessor = 1 << 10,
  NonlocalKernel = 1 << 11,
  NonlocalIntegratedBC = 1 << 12,
  InternalSideIndicator = 1 << 13,
  TransientMultiApp = 1 << 14,
  MultiAppTransfer = 1 << 15
};

template <typename T>
class EnumFlag
{
public:
  using UnderlyingType = typename std::underlying_type<T>::type;
  EnumFlag(const T & flags) : m_flags(static_cast<UnderlyingType>(flags)) {}
  bool operator&(T r) const { return 0 != (m_flags & static_cast<UnderlyingType>(r)); }
  operator Interfaces() const { return Interfaces(m_flags); }
  static const T NoFlag = static_cast<T>(0);

private:
  UnderlyingType m_flags;
};
template <typename T>
EnumFlag<typename std::enable_if<std::is_same<T, Interfaces>::value, T>::type>
operator|(T l, T r)
{
  return static_cast<T>(static_cast<typename EnumFlag<T>::UnderlyingType>(l) |
                        static_cast<typename EnumFlag<T>::UnderlyingType>(r));
}

struct Attribute
{
  AttributeId id;
  int64_t value;
  std::string strvalue;
  inline bool operator==(const Attribute & other) const
  {
    return id == other.id && value == other.value && strvalue == other.strvalue;
  }
  inline bool operator<(const Attribute & other) const
  {
    if (id == other.id)
    {
      if (value == other.value)
        return strvalue < other.strvalue;
      return value < other.value;
    }
    return id < other.id;
  }
};

class TheWarehouse
{
public:
  TheWarehouse();
  ~TheWarehouse();

  class Builder
  {
  public:
    Builder(TheWarehouse & w) : _w(w) {}
    Builder thread(int tid)
    {
      _attribs.push_back({AttributeId::Thread, tid, ""});
      return *this;
    }
    Builder name(const std::string & name)
    {
      _attribs.push_back({AttributeId::Name, 0, name});
      return *this;
    }
    Builder interfaces(int ifaces)
    {
      _attribs.push_back({AttributeId::Interfaces, ifaces, ""});
      return *this;
    }
    Builder interfaces(Interfaces ifaces)
    {
      _attribs.push_back({AttributeId::Interfaces, (int)ifaces, ""});
      return *this;
    }
    Builder subdomain(int id)
    {
      _attribs.push_back({AttributeId::Subdomain, id, ""});
      return *this;
    }
    Builder boundary(int id)
    {
      _attribs.push_back({AttributeId::Boundary, id, ""});
      return *this;
    }
    Builder exec_on(int on)
    {
      _attribs.push_back({AttributeId::ExecOn, on, ""});
      return *this;
    }
    Builder system(const std::string & sys)
    {
      _attribs.push_back({AttributeId::System, 0, sys});
      return *this;
    }

    /// TODO: delete this later - it is a temporary hack for dealing with inter-system dependencies
    Builder pre_ic(bool pre_ic)
    {
      _attribs.push_back({AttributeId::PreIC, (int)pre_ic, ""});
      return *this;
    }
    /// TODO: delete this later - it is a temporary hack for dealing with inter-system dependencies
    Builder pre_aux(bool pre_aux)
    {
      _attribs.push_back({AttributeId::PreAux, (int)pre_aux, ""});
      return *this;
    }
    int prepare() { return _w.prepare(_attribs); }
    size_t count() { return _w.count(prepare()); }
    std::vector<Attribute> attribs() { return _attribs; }
    template <typename T>
    int queryInto(std::vector<T *> & results)
    {
      return _w.queryInto(_attribs, results);
    }

  private:
    TheWarehouse & _w;
    std::vector<Attribute> _attribs;
  };
  Builder build() { return Builder(*this); }

  void add(std::shared_ptr<MooseObject> obj, const std::string & system);

  /// Any attributes specified in extras overwrite/trump ones read from the object's current state.
  void update(MooseObject * obj, const std::vector<Attribute> & extras = {});

  /// prepares a query and returns an associated query_id (i.e. for use with the query function).
  int prepare(const std::vector<Attribute> & conds);

  const std::vector<MooseObject *> query(int query_id);

  size_t count(int query_id);

  template <typename T>
  int queryInto(const std::vector<Attribute> & conds, std::vector<T *> & results)
  {
    int query_id = -1;
    if (_query_cache.count(conds) == 0)
      query_id = prepare(conds);
    else
      query_id = _query_cache[conds];
    queryInto(query_id, results);
    return query_id;
  }

  template <typename T>
  std::vector<T *> & queryInto(int query_id, std::vector<T *> & results)
  {
    std::cout << "    * querying query_id=" << query_id << "\n";
    auto objs = query(query_id);
    results.resize(objs.size());
    for (unsigned int i = 0; i < objs.size(); i++)
    {
      auto obj = objs[i];
      assert(dynamic_cast<T *>(obj));
      results[i] = dynamic_cast<T *>(obj);
    }
    return results;
  }

private:
  void readAttribs(const MooseObject * obj,
                   const std::string & system,
                   std::vector<Attribute> & attribs);

  std::unique_ptr<Storage> _store;
  std::vector<std::shared_ptr<MooseObject>> _objects;
  std::unordered_map<MooseObject *, int> _obj_ids;

  std::vector<std::vector<MooseObject *>> _obj_cache;
  std::map<std::vector<Attribute>, int> _query_cache;
};

#endif // THEWAREHOUSE_H
