[Mesh]
  [./generated]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    nx = 3
    ymin = -2
    ymax = 3
    ny = 5
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./u]
  [../]
[]

[VectorPostprocessors]
  [./face_info]
    type = TestFaceInfo
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
