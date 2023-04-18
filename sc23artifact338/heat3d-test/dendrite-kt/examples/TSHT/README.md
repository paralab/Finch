Running codes (argv and config.txt)
=============
TSHT tutorial consists the following additional features over SSHT:
* Timestepper (Backward Euler and BDF2)
* Internal forcing for 2D and 3D or pure decay.

Running the code 
* from argv (fixed timestepper, no forcing)
```bash
$exec eleOrder level mfree -petsc_config /  $exec 2 4 0 -ksp_type bcgs ...
``` 

* from config.txt
```
basisFunction = "linear" # "quadratic"
mfree = false # "true"
dt = 0.1
totalT = 1.0
BDF2 = "true" # "false"
forcing = "true" # "false"
background_mesh = {
  refine_lvl = 4
  min = [0.0, 0.0]   
  max = [1.0, 1.0]
}

#################### solver setting ####################
solver_options_tsht = {
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "asm"
  ksp_atol = 1e-7
  ksp_rtol = 1e-10
  ksp_converged_reason = ""
  ksp_monitor = ""
}
```

Regression tests
============
Setting from Clion:
```text
* Setting->build->toolchains: create "PETSC" toolchain with C and C++ compiler
* Setting->build->CMake: use the PETSC toolchain and keep the name as "Debug-PETSC-2D" or "Release-PETSC"
```
2D and 3D tests:
* make sure your 2D/3D executable is in one of the following folders (it will choose release over debug if exists):
  * cmake-build-release-petsc-2d / cmake-build-release-petsc
  * cmake-build-release-2d / cmake-build-release
  * cmake-build-debug-petsc-2d / cmake-build-debug-petsc
  * cmake-build-debug-2d / cmake-build-debug

```bash
python3 ssht_regression_test.py --runner=$(which mpirun)
```
Expect small difference (<1e-10) for some np > 1 tests.