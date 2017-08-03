# Using CappedMohrCoulomb with compressive failure only
# checking for small deformation
# A single element is stretched equally in all directions.
# This causes the return direction to be along the sigma_I = sigma_II = sigma_III line
# compressive_strength is set to 1Pa, and smoothing_tol = 0.1Pa
# The smoothed yield function comes from two smoothing operations.
# The first is on sigma_I and sigma_II (sigma_I >= sigma_II >= sigma_III):
# yf = -sigma_I + ismoother(0) - compressive_strength
#    = -sigma_I + (0.5 * smoothing_tol - smoothing_tol / Pi) - compressive_strength
#    = -sigma_I + 0.018169 - 1
# The second has the argument of ismoother equal to -0.018169.
# ismoother(-0.018169) = 0.5 * (-0.018169 + 0.1) - 0.1 * cos (0.5 * Pi * -0.018169 / 0.1) / Pi
#                     = 0.010372
# So the final yield function is
# yf = -sigma_I + 0.018169 + 0.010372 - 1 = -sigma_I + 0.028541 - 1
# However, because of the asymmetry in smoothing (the yield function is obtained
# by first smoothing -sigma_I-cs and -sigma_II-cs, and then by smoothing this
# result with -sigma_III-cs) the result is sigma_I > sigma_II = sigma_III


[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]


[BCs]
  [./x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 'front back'
    function = '-1E-6*x'
  [../]
  [./y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 'front back'
    function = '-1E-6*y'
  [../]
  [./z]
    type = FunctionPresetBC
    variable = disp_z
    boundary = 'front back'
    function = '-1E-6*z'
  [../]
[]

[AuxVariables]
  [./stress_I]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_II]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_III]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./yield_fcn]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_I]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_I
    scalar_type = MaxPrincipal
  [../]
  [./stress_II]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = MidPrincipal
    variable = stress_II
  [../]
  [./stress_III]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = MinPrincipal
    variable = stress_III
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./yield_fcn_auxk]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 0
    variable = yield_fcn
  [../]
[]

[Postprocessors]
  [./s_I]
    type = PointValue
    point = '0 0 0'
    variable = stress_I
  [../]
  [./s_II]
    type = PointValue
    point = '0 0 0'
    variable = stress_II
  [../]
  [./s_III]
    type = PointValue
    point = '0 0 0'
    variable = stress_III
  [../]
  [./s_xx]
    type = PointValue
    point = '0 0 0'
    variable = stress_xx
  [../]
  [./s_xy]
    type = PointValue
    point = '0 0 0'
    variable = stress_xy
  [../]
  [./s_xz]
    type = PointValue
    point = '0 0 0'
    variable = stress_xz
  [../]
  [./s_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./s_yz]
    type = PointValue
    point = '0 0 0'
    variable = stress_yz
  [../]
  [./s_zz]
    type = PointValue
    point = '0 0 0'
    variable = stress_zz
  [../]
  [./f]
    type = PointValue
    point = '0 0 0'
    variable = yield_fcn
  [../]
[]

[UserObjects]
  [./ts]
    type = TensorMechanicsHardeningConstant
    value = 1
  [../]
  [./cs]
    type = TensorMechanicsHardeningConstant
    value = 1
  [../]
  [./coh]
    type = TensorMechanicsHardeningConstant
    value = 1E6
  [../]
  [./ang]
    type = TensorMechanicsHardeningConstant
    value = 0.5
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '0 2.0E6'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./tensile]
    type = CappedMohrCoulombStressUpdate
    tensile_strength = ts
    compressive_strength = cs
    cohesion = coh
    friction_angle = ang
    dilation_angle = ang
    smoothing_tol = 0.1
    yield_function_tol = 1.0E-12
  [../]
  [./stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = tensile
    perform_finite_strain_rotations = false
  [../]
[]


[Executioner]
  end_time = 1
  dt = 1
  type = Transient
[]


[Outputs]
  file_base = small_deform12
  csv = true
[]