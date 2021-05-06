[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./uo_e]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxKernels]
  [./uo_reporter]
    type = MatPropUserObjectAux
    variable = uo_e
    material_user_object = uo2
    execute_on = 'linear timestep_end'
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 1
  [../]

  [./right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 2
  [../]
[]

[Materials]
  [./mat1]
    type = UserObjectMaterial
    mat_prop = 'prop1'
    user_object = uo2
  [../]
  [./mat2]
    type = GenericConstantMaterial
    prop_names = 'prop2'
    prop_values = 2.718281828459
  [../]
[]

[UserObjects]
  [./uo1]
    type = MaterialPropertyUserObject
    mat_prop = prop1
    execute_on = linear
  [../]
  [./uo2]
    type = MaterialPropertyUserObject
    mat_prop = prop2
    execute_on = linear
  [../]
[]

[Executioner]
  type = Steady

  solve_type = 'PJFNK'
[]

[Outputs]
  file_base = uo_material
  exodus = true
[]
