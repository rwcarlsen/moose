/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "SlopeLimitingBase.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include <unistd.h>
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"

template<>
InputParameters validParams<SlopeLimitingBase>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Base class for slope limiting to limit the slopes of cell average variables.");
  params.addParam<bool>("include_bc", true, "Indicate whether to include bc, default = true");
  params.addRequiredParam<UserObjectName>("slope_reconstruction", "Name for slope reconstruction user object");
  return params;
}

SlopeLimitingBase::SlopeLimitingBase(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _include_bc(getParam<bool>("include_bc")),
    _rslope(getUserObject<SlopeReconstructionBase>("slope_reconstruction")),
    _q_point_face(_assembly.qPointsFace()),
    _qrule_face(_assembly.qRuleFace()),
    _JxW_face(_assembly.JxWFace()),
    _normals_face(_assembly.normals()),
    _side(_assembly.side()),
    _side_elem(_assembly.sideElem()),
    _side_volume(_assembly.sideElemVolume()),
    _neighbor_elem(_assembly.neighbor())
{
}

void
SlopeLimitingBase::initialize()
{
  // ElementLoopUserObject::initialize();

  // _lslope.clear();
}

void
SlopeLimitingBase::execute()
{
}

void
SlopeLimitingBase::computeSlopeLimiter()
{
  dof_id_type elementID = _current_elem->id();

  _lslope[elementID] = limitElementSlope();
}

void
SlopeLimitingBase::threadJoin(const UserObject & y)
{
  const SlopeLimitingBase & pps = static_cast<const SlopeLimitingBase &>(y);

  // _lslope.insert(pps._lslope.begin(), pps._lslope.end());
}

const std::vector<RealGradient> &
SlopeLimitingBase::getElementSlope(dof_id_type elementid) const
{
  return _lslope.at(elementid);
}

void
SlopeLimitingBase::serialize(std::string & serialized_buffer)
{
#if 0
  std::ostringstream oss;

  // First store the number of elements to send
  unsigned int size = _interface_elem_ids.size();
  oss.write((char *) &size, sizeof(size));

  for (auto it = _interface_elem_ids.begin(); it != _interface_elem_ids.end(); ++it)
  {
    storeHelper(oss, *it, this);
    storeHelper(oss, _lslope[*it], this);
  }

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
#endif
}

void
SlopeLimitingBase::deserialize(std::vector<std::string> & serialized_buffers)
{
#if 0
  // The input string stream used for deserialization
  std::istringstream iss;

  mooseAssert(serialized_buffers.size() == _app.n_processors(), "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  for (auto rank = decltype(_app.n_processors())(0); rank < serialized_buffers.size(); ++rank)
  {
    if (rank == processor_id())
      continue;

    iss.str(serialized_buffers[rank]);    // populate the stream with a new buffer
    iss.clear();                          // reset the string stream state

    // Load the communicated data into all of the other processors' slots

    unsigned int size = 0;
    iss.read((char *) &size, sizeof(size));

    for (unsigned int i = 0; i < size; i++)
    {
      dof_id_type key;
      loadHelper(iss, key, this);

      std::vector<RealGradient> value;
      loadHelper(iss, value, this);

      // merge the data we received from other procs
      _lslope.insert(std::pair<dof_id_type, std::vector<RealGradient> >(key, value));
    }
  }
#endif
}

void
SlopeLimitingBase::finalize()
{
  // ElementLoopUserObject::finalize();

#if 0
  if (_app.n_processors() > 1)
  {
    std::vector<std::string> send_buffers(1);
    std::vector<std::string> recv_buffers;

    recv_buffers.reserve(_app.n_processors());
    serialize(send_buffers[0]);
    comm().allgather_packed_range((void *)(nullptr), send_buffers.begin(), send_buffers.end(), std::back_inserter(recv_buffers));
    deserialize(recv_buffers);
  }
#endif
}
