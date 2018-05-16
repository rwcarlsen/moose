//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TheWarehouse.h"
#include "MooseObject.h"
#include "TaggingInterface.h"
#include "BoundaryRestrictable.h"
#include "BlockRestrictable.h"
#include "SetupInterface.h"

#include <memory>

class Storage
{
public:
  virtual void add(int obj_id, const std::vector<Attribute> & attribs) = 0;
  virtual std::vector<int> query(const std::vector<Attribute> & conds) = 0;
  virtual void set(int obj_id, const std::vector<Attribute> & attribs) = 0;
};

class VecStore : public Storage
{
private:
  struct Data
  {
    int id;
    std::string type;
    std::string name;
    std::string system;
    int thread = 0;
    bool enabled = true;
    std::vector<boundary_id_type> boundaries;
    std::vector<subdomain_id_type> subdomains;
    std::vector<int> execute_ons;
    std::vector<int> vector_tags;
    std::vector<int> matrix_tags;
  };

public:
  virtual void add(int obj_id, const std::vector<Attribute> & attribs) override
  {
    if (obj_id < _data.size())
      throw std::runtime_error("object with id " + std::to_string(obj_id) + " already added");

    _data.push_back({});
    auto & d = _data.back();
    d.id = obj_id;
    set(d, attribs);
  }

  virtual std::vector<int> query(const std::vector<Attribute> & conds) override
  {
    std::vector<int> objs;
    for (int i = 0; i < _data.size(); i++)
    {
      auto & d = _data[i];
      bool passes = true;
      for (auto & cond : conds)
      {
        switch (cond.id)
        {
          case AttributeId::Thread:
            passes = cond.value == d.thread;
            break;
          case AttributeId::Type:
            passes = cond.strvalue == d.type;
            break;
          case AttributeId::Name:
            passes = cond.strvalue == d.name;
            break;
          case AttributeId::System:
            passes = cond.strvalue == d.system;
            break;
          case AttributeId::Enabled:
            passes = cond.value == d.enabled;
            break;
          case AttributeId::Boundary:
            passes = false;
            for (auto val : d.boundaries)
              if (cond.value == Moose::ANY_BOUNDARY_ID || val == Moose::ANY_BOUNDARY_ID ||
                  cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          case AttributeId::Subdomain:
            passes = false;
            for (auto val : d.subdomains)
              if (cond.value == Moose::ANY_BLOCK_ID || val == Moose::ANY_BLOCK_ID ||
                  cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          case AttributeId::ExecOn:
            passes = false;
            for (auto val : d.execute_ons)
              if (cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          case AttributeId::VectorTag:
            passes = false;
            for (auto val : d.vector_tags)
              if (cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          case AttributeId::MatrixTag:
            passes = false;
            for (auto val : d.matrix_tags)
              if (cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          default:
            throw std::runtime_error("unknown AttributeId " +
                                     std::to_string(static_cast<int>(cond.id)));
        }
        if (!passes)
          break;
      }
      if (passes)
        objs.push_back(i);
    }
    return objs;
  }

  virtual void set(int obj_id, const std::vector<Attribute> & attribs) override
  {
    Data * dat = nullptr;
    for (auto & d : _data)
      if (d.id == obj_id)
      {
        dat = &d;
        break;
      }

    if (!dat)
      throw std::runtime_error("unknown object id " + std::to_string(obj_id));

    set(*dat, attribs);
  }

private:
  void set(Data & d, const std::vector<Attribute> & attribs)
  {
    for (auto & attrib : attribs)
    {
      switch (attrib.id)
      {
        case AttributeId::Thread:
          d.thread = attrib.value;
          break;
        case AttributeId::Type:
          d.type = attrib.value;
          break;
        case AttributeId::Name:
          d.name = attrib.value;
          break;
        case AttributeId::System:
          d.system = attrib.strvalue;
          break;
        case AttributeId::Enabled:
          d.enabled = attrib.value;
          break;
        case AttributeId::Boundary:
          d.boundaries.push_back(attrib.value);
          break;
        case AttributeId::Subdomain:
          d.subdomains.push_back(attrib.value);
          break;
        case AttributeId::ExecOn:
          d.execute_ons.push_back(attrib.value);
          break;
        case AttributeId::VectorTag:
          d.vector_tags.push_back(attrib.value);
          break;
        case AttributeId::MatrixTag:
          d.matrix_tags.push_back(attrib.value);
          break;
        default:
          throw std::runtime_error("unknown AttributeId " +
                                   std::to_string(static_cast<int>(attrib.id)));
      }
    }
  }

  std::vector<Data> _data;
};

TheWarehouse::TheWarehouse() : _store(new VecStore()){};

void
TheWarehouse::add(std::unique_ptr<MooseObject> obj, const std::string & system)
{
  for (int i = 0; i < _query_dirty.size(); i++)
    _query_dirty[i] = true;

  std::vector<Attribute> attribs;
  readAttribs(obj.get(), system, attribs);
  _objects.push_back(std::move(obj));
  int obj_id = _objects.size() - 1;
  _store->add(obj_id, attribs);
}

void
TheWarehouse::update(std::unique_ptr<MooseObject> obj, const std::string & system)
{
  for (int i = 0; i < _query_dirty.size(); i++)
    _query_dirty[i] = true;

  int obj_id = -1;
  for (int i = 0; i < _objects.size(); i++)
    if (_objects[i] == obj)
    {
      obj_id = i;
      break;
    }

  if (obj_id == -1)
    throw std::runtime_error("cannot update unknown object");

  std::vector<Attribute> attribs;
  readAttribs(obj.get(), system, attribs);
  _store->set(obj_id, attribs);
}

// prepares a query and returns an associated query_id (i.e. for use with the query function.
int
TheWarehouse::prepare(const std::vector<Attribute> & conds)
{
  auto obj_ids = _store->query(conds);
  _query_dirty.push_back(true);
  _obj_cache.push_back({});
  _query_cache.push_back(conds);
  return _obj_cache.size() - 1;
}

const std::vector<MooseObject *> &
TheWarehouse::query(int query_id)
{
  if (query_id >= _obj_cache.size())
    throw std::runtime_error("unknown query id");

  if (_query_dirty[query_id])
  {
    auto & vec = _obj_cache[query_id];
    vec.clear();
    for (auto & id : _store->query(_query_cache[query_id]))
      vec.push_back(_objects[id].get());
    _query_dirty[query_id] = false;
  }

  return _obj_cache[query_id];
}

void
TheWarehouse::readAttribs(MooseObject * obj,
                          const std::string & system,
                          std::vector<Attribute> & attribs)
{
  attribs.push_back({AttributeId::System, 0, system});
  attribs.push_back({AttributeId::Type, 0, ""}); // TODO: make this not be blank
  attribs.push_back({AttributeId::Name, 0, obj->name()});
  attribs.push_back({AttributeId::Thread, static_cast<int>(obj->getParam<THREAD_ID>("_tid")), ""});
  attribs.push_back({AttributeId::Enabled, obj->enabled(), ""});

  auto ti = dynamic_cast<TaggingInterface *>(obj);
  if (ti)
  {
    for (auto tag : ti->getVectorTags())
      attribs.push_back({AttributeId::VectorTag, static_cast<int>(tag), ""});
    for (auto tag : ti->getMatrixTags())
      attribs.push_back({AttributeId::MatrixTag, static_cast<int>(tag), ""});
  }
  auto blk = dynamic_cast<BlockRestrictable *>(obj);
  if (blk)
  {
    for (auto id : blk->blockIDs())
      attribs.push_back({AttributeId::Subdomain, id, ""});
  }
  auto bnd = dynamic_cast<BoundaryRestrictable *>(obj);
  if (bnd)
  {
    for (auto & bound : bnd->boundaryIDs())
      attribs.push_back({AttributeId::Boundary, bound, ""});
  }
  auto sup = dynamic_cast<SetupInterface *>(obj);
  if (sup)
  {
    for (auto & on : sup->getExecuteOnEnum().items())
      attribs.push_back({AttributeId::ExecOn, on, ""});
  }
}
