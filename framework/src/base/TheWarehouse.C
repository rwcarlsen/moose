//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

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
    std::string id;
    std::string type;
    std::string name;
    std::string system;
    int thread = 0;
    bool enabled = true;
    std::vector<int> boundaries;
    std::vector<int> subdomains;
    std::vector<int> execute_ons;
    std::vector<std::string> tags;
  };

public:
  virtual void add(int obj_id, const std::vector<Attribute> & attribs) override;
  {
    if (obj_id < _system.size())
      throw std::runtime_error("object with id " + std::to_string(obj_id) + " already added");

    _data.push_back({});
    auto & d = _data.back();
    d.id = obj_id;
    set(d, attribs);
  }

  virtual std::vector<int> query(const std::vector<Attribute> & conds) override;
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
          case AttributeId::System:
            passes = cond.strvalue == d.system;
            break;
          case AttributeId::Enabled:
            passes = cond.value == d.enabled;
            break;
          case AttributeId::Boundary:
            passes = false;
            for (auto val : d.boundaries)
              if (cond.value == val)
              {
                passes = true;
                break;
              }
            break;
          case AttributeId::Subdomain:
            passes = false;
            for (auto val : d.subdomains)
              if (cond.value == val)
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
          case AttributeId::Tag:
            passes = false;
            for (auto val : d.tags)
              if (cond.strvalue == val)
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

  virtual void set(int obj_id, const std::vector<Attribute> & attribs) override;
  {
    Data * d = nullptr;
    for (auto & d : _data)
      if (d.id == obj_id)
      {
        d = &d;
        break;
      }

    if (!d)
      throw std::runtime_error("unknown object id " std::to_string(obj_id));

    set(*d, attribs);
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
        case AttributeId::Tag:
          d.tags.push_back(attrib.strvalue);
          break;
        default:
          throw std::runtime_error("unknown AttributeId " +
                                   std::to_string(static_cast<int>(attrib.id)));
      }
    }
  }

  std::vector<Data> _data;
};

TheWarehouse::TheWarehouse() : _store(nullptr) { _store = libmesh_make_unique<VecStore>(); };

void
TheWarehouse::add(std::unique_ptr<MooseObject> obj)
{
  for (int i = 0; i < _query_dirty.size(); i++)
    _query_dirty[i] = true;

  std::vector<Attribute> attribs;
  readAttribs(attribs);
  _objects.push_back(std::move(obj));
  int obj_id = _objects.size() - 1;
  obj->parameters().set<int>("_warehouse_id") = obj_id;
  _store->add(obj_id, attribs);
}

void
TheWarehouse::readAttribs(MooseObject * obj, std::vector<Attribute> & attribs)
{
  attribs.push_back({AttributeId::System, 0, obj->system});
  attribs.push_back({AttributeId::Thread, obj->thread, ""});
  attribs.push_back({AttributeId::Enabled, obj->enabled, ""});
  for (auto & tag : obj->tags)
    attribs.push_back({AttributeId::Tag, 0, tag});
  for (auto & sub : obj->subdomains)
    attribs.push_back({AttributeId::Subdomain, sub, ""});
  for (auto & bound : obj->boundaries)
    attribs.push_back({AttributeId::Boundary, bound, ""});
  for (auto & on : obj->execute_ons)
    attribs.push_back({AttributeId::ExecOn, on, ""});
}

TheWarehouse::update(std::unique_ptr<MooseObject> obj)
{
  for (int i = 0; i < _query_dirty.size(); i++)
    _query_dirty[i] = true;

  int obj_id = obj->getParam<int>("_warehouse_id");
  std::vector<Attribute> attribs;
  readAttribs(attribs);
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
