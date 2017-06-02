
class ResidualBuffer
{
public:
  ResidualBuffer(DofMap & dof_map, NumericVector<Number> & residual_vector)
    : _max_size(0), _residual(residual_vector), _dof_map(dof_map)
  {
  }

  void flush()
  {
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _buffer.size(); i++)
        _residual.add(_buffer[i].row, _buffer[i].val);
    }
    if (_max_size < _buffer.size())
      _max_size = _buffer.size();

    _buffer.clear();
    _buffer.reserve(_max_size * 2);
  }

  ~ResidualBuffer() { flush(); }

  void customContribution(dof_id_type dof, Real value) { _buffer.push_back({dof, value}); }

  void varContribution(DenseVector<Number> & res_block,
                       std::vector<dof_id_type> & dof_indices,
                       Real scaling_factor)
  {
    if (dof_indices.size() == 0 || res_block.size() == 0)
      return;

    _dof_map.constrain_element_vector(res_block, dof_indices, false);

    if (scaling_factor != 1.0)
      for (unsigned int i = 0; i < res_block.size(); i++)
        _buffer.push_back({_temp_dof_indices[i], res_block[i]});
    else
      for (unsigned int i = 0; i < res_block.size(); i++)
        _buffer.push_back({_temp_dof_indices[i], res_block[i] * scaling_factor});
  }

private:
  struct Entry
  {
    dof_id_type row;
    Real val;
  }

  unsigned int _max_size;
  NumericVector<Number> & _residual;
  DofMap & _dof_map;
  std::vector<Entry> _buffer;
};

class AutoBuffer
{
public:
  AutoBuffer(FEProblemBase & prob, THREAD_ID tid)
    : _sys(prob.getNonlinearSystem()), _assem(prob.assembly()), _tid(tid)
  {
    for (auto type : std::vector({KT_TIME, KT_NOTIME}))
      if (_sys.hasResidualVector(type))
      {
        _types.push_back(type);
        _buffs.push_back(ResidualBuffer(_sys.dofMap(), _sys.residualVector(type)));
      }
  }

  void recordResiduals()
  {
    const std::vector<MooseVariable *> & vars = _sys.getVariables(_tid);
    for (const auto & var : vars)
      for (unsigned int i = 0; i < _buffs.size(); i++)
      {
        auto & block = _assem.residualBlock(var, _types[i]);
        _buffs[i].varContribution(block, var->dof_indices(), var->scalingFactor());
      }
  }

  void flush()
  {
    for (auto & buf : _buffs)
      buf.flush();
  }

private:
  SystemBase & _sys;
  Assembly & _assem;
  THREAD_ID _tid;
  std::vector<Moose::KernelType> _types;
  std::vector<ResidualBuffer> _buffs;
};
