#!/usr/bin/env python3

import os
import sys
import re

this_dir = os.path.dirname(os.path.realpath(__file__))
scripts_dir = os.path.join(this_dir, '..', '..', 'scripts')
sys.path.append(scripts_dir)

# choose release mode when possible
found_2d_exec, found_3d_exec = False, False
for build_folder in ['cmake-build-debug-petsc', 'cmake-build-debug']:#, 'cmake-build-release-petsc', 'cmake-build-release']:
  exec_2d = os.path.join(this_dir, '..', '..', build_folder + "-2d", 'examples', 'SSHT', 'steady-heat')
  exec_3d = os.path.join(this_dir, '..', '..', build_folder, 'examples', 'SSHT', 'steady-heat')
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

base_cfg = {}

cases = [
    # check serial and parallel
    {'dim': 3, 'order': "linear", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 1},
    {'dim': 3, 'order': "linear", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 4},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 1},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 4},
    {'dim': 3, 'order': "linear", 'refine_lvl': 3, 'mfree': "true", 'ntasks': 1},
    {'dim': 3, 'order': "linear", 'refine_lvl': 3, 'mfree': "true", 'ntasks': 4},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 3, 'mfree': "true", 'ntasks': 1},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 3, 'mfree': "true", 'ntasks': 4},
    #'dim': 3,  check the slope of linear basis functions
    {'dim': 3, 'order': "linear", 'refine_lvl': 2, 'mfree': "false", 'ntasks': 2},
    {'dim': 3, 'order': "linear", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 4},
    {'dim': 3, 'order': "linear", 'refine_lvl': 4, 'mfree': "false", 'ntasks': 4},
    {'dim': 3, 'order': "linear", 'refine_lvl': 5, 'mfree': "false", 'ntasks': 8},
    #'dim': 3,  check the slope of quadratic basis functions
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 2, 'mfree': "false", 'ntasks': 2},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 3, 'mfree': "false", 'ntasks': 4},
    {'dim': 3, 'order': "quadratic", 'refine_lvl': 4, 'mfree': "false", 'ntasks': 4},
    # {'dim': 3, 'order': "quadratic", 'refine_lvl': 5, 'mfree': "false", 'ntasks': 8},
]
for i in range(len(cases)):
  if cases[i]["mfree"] == "true":
    cases[i]["pc"] = "none"
  else:
    cases[i]["pc"] = "asm"

from regression import RegressionTester, RegexDiffMetric, VecDiffMetric

reg = RegressionTester()

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
reg.add_metric(RegexDiffMetric('Final L2 error', 'output.txt', re.compile(r'L2Error : (.*)'), 5e-8))
reg.add_file_generator("config.txt", """
basisFunction = "{order}"
mfree = {mfree}
background_mesh = {{
  refine_lvl = {refine_lvl}
  min = [0, 0, 0]
  max = [1, 1, 1]
}}

dump_vec = true

#################### solver setting ####################
solver_options_ssht = {{
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "{pc}"
  ksp_atol = 1e-7
  ksp_rtol = 1e-10
  ksp_converged_reason = ""
  ksp_monitor = ""
}}
""")

reg.set_run_cmd(executable_3d)
reg.run_main()

# 2D tests
reg = None
reg = RegressionTester()

for i in range(len(cases)):
  cases[i]["dim"] = 2
  if cases[i]["mfree"] == "true":
    cases[i]["pc"] = "none"
  else:
    cases[i]["pc"] = "asm"

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
reg.add_metric(RegexDiffMetric('Final L2 error', 'output.txt', re.compile(r'L2Error : (.*)'), 5e-8))
reg.add_file_generator("config.txt", """
basisFunction = "{order}"
mfree = {mfree}
background_mesh = {{
  refine_lvl = {refine_lvl}
  min = [0, 0]
  max = [1, 1]
}}

dump_vec = true

#################### solver setting ####################
solver_options_ssht = {{
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "{pc}"
  ksp_atol = 1e-7
  ksp_rtol = 1e-10
  ksp_converged_reason = ""
  ksp_monitor = ""
}}
""")
reg.set_run_cmd(executable_2d)
reg.run_main()
