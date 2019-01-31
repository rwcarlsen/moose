
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
calcPropertyDerivs(DualReal & n, Real val, Real dudq)
{
  n.value() = val;
  for (auto& deriv : n.derivatives())
    deriv *= dudq;
}

void
calcPropertyDerivs(Real & n, Real val, Real /*dudq*/)
{
  n = val;
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
  ppsat
};

// only these variables can be involved in identifying a state. And we only need partial
// derivatives w.r.t. the first variable in each pair.
enum class Combo
{
  // imply pT formulation
  pT,
  ph,
  ps,
  pv,
  // imply ve formulation
  vp,
  ve,
  vh,
  vT,
};

// clang-format off
template <ComputeStage compute_stage>
class SinglePhaseState
{
  SinglePhaseState(Combo given, ADReal v1, ADReal v2)
    : _given(given), _prop1(v1), _prop2(v2)
  {
    switch (given)
    {
    case pT:
    case ph:
    case ps:
    case pv:
      _given1 = Var::p;
      break;
    case vp:
    case ve:
    case vh:
    case vT:
      _given1 = Var::v;
      break;
    }
  }

  #define propfunc(propvar) \
  ADReal propvar() \
  { \
    ADReal v = _prop1; \
    calcPropertyDeriv(v, _##propvar(), d##propvar(_given1)); \
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

  Real prop1() {return _prop1.value();}
  Real prop2() {return _prop2.value();}

private:
  Real _rho() { return 1 / _v(); } // special case

  ADReal _prop1;
  ADReal _prop2;
  Combo _given;
  Var _given1;
};
// clang-format on

template <ComputeStage compute_stage>
class LegacyWrapper : public SinglePhaseState<compute_stage>
{
public:
  LegacyWrapper(SinglePhaseFluidProperties & spp, Combo g, ADReal v1, ADReal v2)
    : SinglePhaseState<compute_stage>(g, v1, v2), _spp(spp)
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
        return 1;
      case Var::v:
        _spp.p_from_v_e(_v(), _e(), dummy1, dudq, dummy2);
        return dudq;
      default:
        mooseError("fluidprops state does not support p derivative w.r.t. var ", dvar);
    }
  }
  Real dT(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dummy3 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.T_from_v_e(_v(), _e(), dummy1, dudq, dummy2);
        _spp.v_from_p_T(_p(), _T(), dummy1, dummy2, dummy3);
        dudq *= dummy2;
        return dudq;
      case Var::v:
        _spp.T_from_v_e(_v(), _e(), dummy1, dudq, dummy2);
        return dudq;
      case Var::T:
        return 1;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real de(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.e_from_p_rho(_p(), _rho(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        _spp.e_from_v_h(_v(), _h(), dummy1, dudq, dummy2);
        return dudq;
      case Var::e:
        return 1;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real drho(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.rho_from_p_T(_p(), _T(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        return -1/_v()/_v();
      case Var::rho:
        return 1;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real dh(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.h_from_p_T(_p(), _T(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        _spp.h_from_T_v(_T(), _v(), dummy1, dummy1, dudq);
        return dudq;
      case Var::h:
        return 1;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real ds(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.s_from_p_T(_p(), _T(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        _spp.s_from_v_e(_T(), _v(), dummy1, dudq, dummy2);
        return dudq;
      case Var::s:
        return 1;
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real dv(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      case Var::p:
        _spp.v_from_p_T(_p(), _T(), dummy1, dudq, dummy2);
        return dudq;
      case Var::v:
        return 1;
      case Var::rho:
        return -1/_rho()/_rho();
      default:
        mooseError("fluidprops state does not support T derivative w.r.t. var ", dvar);
    }
  }
  Real dcp(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      default:
        mooseError("fluidprops state does not support cp derivative w.r.t. var ", dvar);
    }
  }
  Real dcv(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      default:
        mooseError("fluidprops state does not support cv derivative w.r.t. var ", dvar);
    }
  }
  Real dmu(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      default:
        mooseError("fluidprops state does not support mu derivative w.r.t. var ", dvar);
    }
  }
  Real dg(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      default:
        mooseError("fluidprops state does not support g derivative w.r.t. var ", dvar);
    }
  }
  Real dbeta(Var dvar)
  {
    Real dummy1 = 0, dummy2 = 0, dudq = 0;
    switch (dvar)
    {
      default:
        mooseError("fluidprops state does not support beta derivative w.r.t. var ", dvar);
    }
  }

  // clang-format off
  Real _p()
  {
    switch (_given)
    {
      case Combo::pT:
      case Combo::ph:
      case Combo::ps:
      case Combo::pv:
        return prop1();
      case Combo::vp:
        return prop2();
      case ve:
      case vh:
      case vT:
        return _spp.p_from_v_e(prop1(), _e());
    }
  }
  Real _e()
  {
    switch (_given)
    {
      case Combo::ve:
        return prop2();
      case Combo::ps:
      case Combo::pT:
      case Combo::ph:
      case Combo::pv:
        return _spp.e_from_p_T(prop1(), _T());
      case Combo::vp:
      case Combo::vh:
      case Combo::vT:
        return _spp.e_from_T_v(_T(), prop1());
    }
  }
  Real _h()
  {
    switch (_given)
    {
      case Combo::ph:
      case Combo::vh:
        return prop2();
      case Combo::pT:
      case Combo::ps:
      case Combo::pv:
        return _spp.h_from_p_T(prop1(), _T());
      case Combo::vp:
      case Combo::ve:
      case Combo::vT:
        return _spp.h_from_p_T(_p(), _T());
    }
  }
  Real _T()
  {
    switch (_given)
    {
      case Combo::vT:
      case Combo::pT:
        return prop2();
      case Combo::ph:
      case Combo::ps:
      case Combo::pv:
        return _spp.T_from_v_e(_v(), _e());
      case Combo::vh:
      case Combo::vp:
      case Combo::ve:
        return _spp.T_from_v_e(prop1(), _e());
    }
  }
  Real _s()
  {
    switch (_given)
    {
      case Combo::ps:
        return prop2();
      case Combo::pT:
      case Combo::ph:
      case Combo::pv:
        return _spp.s_from_p_T(prop1(), _T());
      case Combo::vp:
      case Combo::ve:
      case Combo::vh:
      case Combo::vT:
        return _spp.s_from_v_e(prop1(), _e());
    }
  }
  Real _v()
  {
    switch (_given)
    {
      case Combo::pT:
      case Combo::ph:
      case Combo::ps:
        return _spp.v_from_p_T(prop1(), _T());
      case Combo::pv:
        return prop2();
      case Combo::vT:
      case Combo::vh:
      case Combo::vp:
      case Combo::ve:
        return prop1();
    }
  }
  Real _cp();
  Real _cv();
  Real _mu();
  Real _g();
  Real _beta();
#undef propcase
  // clang-format on
private:
  const SinglePhaseFluidProperties & _spp;
};
} // namespace fluidprops

