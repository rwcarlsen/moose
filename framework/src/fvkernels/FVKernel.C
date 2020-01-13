
#include "FVKernel.h"

FVKernel::FVKernel(const InputParameters & params) : MooseObject(params), TaggingInterface(this) {}

FVKernelFace::FVKernelFace(const InputParameters & params)
  : FVKernel(params), _matprop_iface(this, {}, {})
{
}
