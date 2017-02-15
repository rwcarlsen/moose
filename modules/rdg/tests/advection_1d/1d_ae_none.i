[GlobalParams]
  order = CONSTANT
  family = MONOMIAL
  slope_reconstruction = rslope
  slope_limiting = lslope
  implicit = false
[]

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 100000
[]

[Problem]
  type = RDGProblem
[]

[Functions]
  [./ic_u]
    type = PiecewiseConstant
    axis = 0
    direction = right
    xy_data = '0.1 0.5
               0.6 1.0
               1.0 0.5'
  [../]
[]

[UserObjects]
  [./rslope]
    type = AEFVSlopeReconstructionOneD
  [../]

  [./lslope]
    type = AEFVSlopeLimitingOneD
    u = u
    scheme = 'none' #none | minmod | mc | superbee
  [../]

  [./internal_side_flux]
    type = AEFVUpwindInternalSideFlux
    u = u
  [../]

  [./free_outflow_bc]
    type = AEFVFreeOutflowBoundaryFlux
    u = u
    boundary = 'left right'
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    variable = 'u'
    function = ic_u
  [../]
[]

[Kernels]
  [./time_u]
    implicit = true
    type = TimeDerivative
    variable = u
  [../]
[]

[DGKernels]
  [./concentration]
    type = RDGFlux
    variable = u
    component = 0
    flux = internal_side_flux
  [../]
[]

[BCs]
  [./concentration]
    type = RDGFluxBC
    boundary = 'left right'
    variable = u
    component = 0
    flux = free_outflow_bc
  [../]
[]

[Materials]
  [./aefv]
    type = AEFVMaterial
    block = 0
    u = u
  [../]
[]

[Executioner]
  type = Transient
  #[./TimeIntegrator]
  #  type = ExplicitMidpoint
  #[../]
  #solve_type = 'LINEAR'
  scheme = explicit-euler

  l_tol = 1e-4
  nl_rel_tol = 1e-20
  nl_abs_tol = 1e-8
  nl_max_its = 60

  start_time = 0.0
  num_steps = 1
  dt = 5e-4
  dtmin = 1e-4
[]

[Outputs]
  [./exodus]
    type = Exodus
      file_base = 1d_ae_none
  [../]
  print_perf_log = true
[]
