#pragma once

class FVFluxKernel : public MooseObject,
                     public TaggingInterface,
                     public TwoMaterialPropertyInterface,
                     public NeighborCoupleable
{
public:
  FVFluxKernel(const InputParameters & params);
  ADReal computeResidual()
  {
    auto r = _face_area * computeQpResidual(type);

    if (ownElement())
    {
      prepareVectorTag(_assembly, _var.number());
      _local_re(0) = r;
      accumulateTaggedLocalResidual();
    }
    if (ownNeighbor())
    {
      prepareVectorTagNeighbor(_assembly, _var.number());
      _local_re(0) = -r;
      accumulateTaggedLocalResidual();
    }
  }
};

