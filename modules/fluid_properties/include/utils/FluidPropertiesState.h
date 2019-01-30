
// Example code for calculating ADReal for fluid props:
//
//    template <ComputeStage compute_stage>
//    ADReal
//    MyFluidPropClass::p_from_v_t(Real v, Real t, unsigned int qp)
//    {
//      ADReal p = ...;
//      auto dpdv = ...;
//      auto dpdt = ...;
//      calcPropertyDerivs(qp, p, dpdv, _volume_var);
//      calcPropertyDerivs(qp, p, dpdt, _temperature_var);
//      return p;
//    }
//
// where _volume_var and _temperature_var must have been appropriately specified/set earlier by
// the user of this class.
//
// This function needs to go in some generic utils location and be defined in a header for many
// objects to use.
//
void
calcPropertyDerivs(unsigned int qp, DualReal & n, Real dudq, const MooseVariable & var)
{
  auto offset = var.number() * var.sys().getMaxVarNDofsPerElem();
  auto & phi = var.phi();
  for (size_t i = 0; i < var.numberOfDofs(); i++)
  {
    auto phi_i = phi[i][qp];
    n.derivatives()[offset + i] = phi_i * dudq;
  }
}

void
calcPropertyDerivs(unsigned int /*qp*/, Real & /*n*/, Real /*dudq*/, const MooseVariable & /*var*/)
{
}

namespace fluidprops
{

enum class Var
{
  p,
  T,
  v,
  e,
  rho,
  h,
  s,
  cp,
  cv,
  beta,
  g,
  mu
};
enum class Combo
{
  pT,
  ph,
  ps,
  prho,
  ve,
  vh,
  Tv,
  hs
};

// clang-format off
template <ComputeStage compute_stage>
class SinglePhaseState
{
public:
  void setPrimalVar(Var prop, MooseVariable & var) { _vars[prop] = var;}

  #define propfunc(propvar) \
  ADReal propvar(unsigned int qp) \
  { \
    ADReal v = _##propvar(); \
    for (auto it = _vars.begin(); it != _vars.end(); ++it) \
      if (it.first == Var::rho) \
        calcPropertyDeriv(qp, v, d##propvar(Var::v)*(-1/_rho()/_rho()), it.second); \
      else \
        calcPropertyDeriv(qp, v, d##propvar(it.first), it.second); \
    return v; \
  }
  propfunc(p)
  propfunc(T)
  propfunc(e)
  propfunc(rho)
  propfunc(h)
  propfunc(s)
  propfunc(v)
  propfunc(cp)
  propfunc(cv)
  propfunc(mu)
  propfunc(g)
  propfunc(beta)
  #undef propfunc

protected:
  // don't need to implement case where dvar==Var::rho
  virtual Real dp(Var dvar) { mooseError("not implemented"); }
  virtual Real dT(Var dvar) { mooseError("not implemented"); }
  virtual Real de(Var dvar) { mooseError("not implemented"); }
  virtual Real drho(Var dvar) { mooseError("not implemented"); }
  virtual Real dh(Var dvar) { mooseError("not implemented"); }
  virtual Real ds(Var dvar) { mooseError("not implemented"); }
  virtual Real dv(Var dvar) { mooseError("not implemented"); }

  virtual Real _p() { mooseError("not implemented"); }
  virtual Real _e() { mooseError("not implemented"); }
  virtual Real _h() { mooseError("not implemented"); }
  virtual Real _s() { mooseError("not implemented"); }
  virtual Real _v() { mooseError("not implemented"); }
  virtual Real _cp() { mooseError("not implemented"); }
  virtual Real _cv() { mooseError("not implemented"); }
  virtual Real _mu() { mooseError("not implemented"); }
  virtual Real _g() { mooseError("not implemented"); }
  virtual Real _beta() { mooseError("not implemented"); }

  /// Universal gas constant (J/mol/K)
  const Real _R = 8.3144598;
  /// Conversion of temperature from Celsius to Kelvin
  const Real _T_c2k = 273.15;

private:
  Real _rho() { return 1 / _v(); } // special case
  std::unordered_map<Var, MooseVariable *> _vars;
};
// clang-format on

template <ComputeStage compute_stage>
class LegacyWrapper : public SinglePhaseState<compute_stage>
{
public:
  LegacyWrapper(SinglePhaseFluidProperties & spp, Combo g, Real v1, Real v2)
    : _spp(spp), _given(g), _prop1(v1), _prop2(v2)
  {
  }

protected:
  Real dp(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    // p, T, v, e, rho, h, s
    switch (dvar)
    {
      case Var::p:
        dudq = 1;
        return dudq;
      case Var::T:
        _spp.p_from_T_v(T(), v(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        _spp.p_from_v_e(v(), e(), dummy1, dudq, dummy2);
        return dudq;
      case Var::e:
        _spp.p_from_v_e(v(), e(), dummy1, dummy2, dudq);
        return dudq;
      case Var::s:
        _spp.p_from_h_s(h(), s(), dummy1, dummy2, dudq);
        return dudq;
      default:
        mooseError("fluidprops state does not support p derivative w.r.t. var ", dvar);
    }
  }
  Real dT(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        return dudq;
      case Var::T:
        return dudq;
      case Var::v:
        return dudq;
      case Var::e:
        return dudq;
      case Var::h:
        return dudq;
      case Var::s:
        return dudq;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real de(Var dvar);
  Real drho(Var dvar);
  Real dh(Var dvar);
  Real ds(Var dvar);
  Real dv(Var dvar);
  Real dcp(Var dvar);
  Real dcv(Var dvar);
  Real dmu(Var dvar);
  Real dg(Var dvar);
  Real dbeta(Var dvar);

#define propcase(want, from1, from2)                                                               \
  case Combo::from1##from2:                                                                        \
    return _spp.want##_from_##from1##_##from2(_prop1, _prop2)
  Real _p()
  {
    switch (_given)
    {
      case Combo::pT:
      case Combo::ph:
      case Combo::ps:
      case Combo::prho:
        return _prop1;
        propcase(p, v, e);
        propcase(p, v, h);
        propcase(p, T, v);
        propcase(p, h, s);
    }
  }
  Real _e();
  Real _h();
  Real _s();
  Real _v();
  Real _cp();
  Real _cv();
  Real _mu();
  Real _g();
  Real _beta();
#undef propcase

private:
  const SinglePhaseFluidProperties & _spp;
  Real _prop1;
  Real _prop2;
  Combo _given;
};
} // namespace fluidprops

