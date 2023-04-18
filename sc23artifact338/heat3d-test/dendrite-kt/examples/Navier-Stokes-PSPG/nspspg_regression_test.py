#!/usr/bin/env python3

import os
import sys
import re

this_dir = os.path.dirname(os.path.realpath(__file__))
scripts_dir = os.path.join(this_dir, '..', '..', 'scripts')
sys.path.append(scripts_dir)

# choose release mode when possible
found_2d_exec, found_3d_exec = False, False
for build_folder in ['cmake-build-debug-petsc', 'cmake-build-debug', 'cmake-build-release-petsc', 'cmake-build-release']:
  exec_2d = os.path.join(this_dir, '..', '..', build_folder + "-2d", 'examples', 'Navier-Stokes-PSPG', 'ns-pspg')
  exec_3d = os.path.join(this_dir, '..', '..', build_folder, 'examples', 'Navier-Stokes-PSPG', 'ns-pspg')
  if os.path.exists(exec_2d):
    found_2d_exec = True
    executable_2d = exec_2d
  if os.path.exists(exec_3d):
    found_3d_exec = True
    executable_3d = exec_3d

assert found_2d_exec, "Program executable missing ({})".format(executable_2d)
assert found_3d_exec, "Program executable missing ({})".format(executable_3d)
print("Executable_2d: {}".format(executable_2d))
print("Executable_3d: {}".format(executable_3d))

base_cfg = {"dim": 3, "V": "", "mfree": "false", "bcl": 0, "totalT": 1, "Re": 100, "MMS": "true", "Formulation": "PSPG"}

config_template = """
################### mesh setting ######################
basisFunction = "{order}"
mfree = {mfree}
background_mesh = {{
  baseLevel = {basel}
  refineLevelBoundary = {bcl}
  min = [0, 0, 0]
  max = [1, 1, 1]
}}

################### timestepper setting ###############
dt{V} = {dt}
totalT{V} = {totalT}
TimeStepper = "{TS}";

################### solver setting ####################
Formulation = "{Formulation}";
Re = {Re};
MMS = {MMS};

dump_vec = true
################### petsc setting #####################
solver_options = {{
  snes_atol = 1e-12
  snes_rtol = 1e-12
  snes_stol = 1e-10
  snes_max_it = 40
  snes_max_funcs = 80000
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "asm"
  # monitor
  snes_monitor = ""
  snes_converged_reason = ""
}};
"""

cases = [
    # check LDC
    {'order': "quadratic", 'basel': 2, 'bcl': 3, 'V': '_V', 'dt': "[1, 0.1]", 'totalT': "[2, 3]", 'TS': 'BDF2', "MMS": "false",'ntasks': 8},
]



from regression import RegressionTester, RegexDiffMetric, VecDiffMetric

reg = RegressionTester()

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
# reg.add_metric(RegexDiffMetric('L2 error each step', 'output.txt', re.compile(r'L2Error : (.*) (.*) (.*) '), 5e-8))
reg.add_file_generator("config.txt", config_template)

reg.set_run_cmd(executable_3d)
reg.run_main()

# 2D tests
reg = None
reg = RegressionTester()
cases = [
    # check MMS, mfree, timestepper, basis function, mpi
    {'order': "linear", 'basel': 3, 'dt': 0.5, 'TS': "BE", 'ntasks': 1},
    {'order': "linear", 'basel': 3, 'dt': 0.5, 'TS': "BE", 'ntasks': 4},
    {'order': "linear", 'basel': 3, 'dt': 0.5, 'TS': "BE", 'mfree': 'true', 'ntasks': 4},
    {'order': "quadratic", 'basel': 3, 'dt': 0.5, 'TS': "BE", 'ntasks': 1},
    {'order': "quadratic", 'basel': 3, 'dt': 0.5, 'TS': "BE", 'ntasks': 4},
    {'order': "linear", 'basel': 5, 'dt': 0.2, 'TS': "BE", 'ntasks': 8},
    {'order': "linear", 'basel': 5, 'dt': 0.2, 'TS': "BDF2", 'ntasks': 8},
    {'order': "linear", 'basel': 5, 'dt': 0.2, 'TS': "CN", 'ntasks': 8},
    # check LDC, check refine,
    {'order': "quadratic", 'basel': 3, 'bcl': 5, 'dt': 10, 'totalT': 30, 'TS': 'BDF2', "MMS": "false",'ntasks': 8},
    # didn't cover: formulation, Re
]
for i in range(len(cases)):
  cases[i]["dim"] = 2

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
reg.add_metric(RegexDiffMetric('L2 error each step', 'output.txt', re.compile(r'L2Error : (.*) (.*) (.*) '), 5e-8))
reg.add_file_generator("config.txt", config_template)
reg.set_run_cmd(executable_2d)
reg.run_main()
