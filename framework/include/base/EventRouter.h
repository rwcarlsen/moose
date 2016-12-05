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

#ifndef EVENTROUTER_H
#define EVENTROUTER_H

class EventRouter
{
public:
  typedef unsigned int EventId;
  typedef unsigned int Id;
  EventRouter(InputParameters& params) : _id(0), _events() {}

  // reservations, registrations, etc. don't need to explicitly support
  // restart/recovery as long as users are performing event-related activities
  // (other than triggering is not limited to ) in their constructors other
  // than triggering - which should generally NOT be done in constructors.

  // reserve generates a unique id for the given event name that must be used
  // when triggering this event.  This is an optional safety mechanism to
  // prevent triggering of events from other than the original, intended
  // location.  Calls to reserve should generally be done in constructors in
  // order to prevent unintended duplicate reservations (which would cause
  // an error).
  EventId reserve(std::string event)
  {
    if (_event_id.count(event) != 0)
      mooseError("event '" + event + "' has already been reserved");
    auto id = _next_event_id;
    _next_event_id++;
    _event_id[event] = id;
    _id_event[id] = event;
    return id;
  }

  // trigger invokes every callback function registered for the given event
  // id.  This should generally NOT be called in constructors.
  void trigger(EventId id) {
    auto event = _id_event[id];
    for (auto& kv : _events[event])
      kv.second(event);
  }


  // trigger invokes every callback function registered to the given event
  // name.  This causes an error if the event name has previously been
  // reserved (in which case you must trigger with the associated EventId).
  // This should generally NOT be called in constructors.
  void trigger(std::string event)
  {
    if (_event_id.count(event) != 0)
      mooseError("event '" + event + "' has already been reserved");
    for (auto& kv : _events[event])
      kv.second(event);
  }

  // registerCallback ensures func is called every time event is triggered.
  // An ID is returned that can be used to unregister for the callbacks.
  // Registrations should generally be done in constructors in order to avoid
  // unwanted duplicate registrations.
  Id registerCallback(std::string event, std::function<(std::string)> func)
  {
    auto id = _id;
    _id++;
    _events[event][id] = func;
    return id;
  }

  // unregister removes the previously registered callback for the given event
  // name and callback id.
  void unregister(std::string event, Id id) { _events[event].erase(id); }

private:
  Id _next_func_id;
  EventId _next_event_id;
  std::map<std::string, EventId> _event_id;
  std::map<EventId, std::string> _id_event;
  std::map<std::string, std::map<Id, std::function<(std::string)>>> _events;
};

#endif /* EVENTROUTER_H */

