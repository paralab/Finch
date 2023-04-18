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
  exec_2d = os.path.join(this_dir, '..', '..', build_folder + "-2d", 'examples', 'TSHT', 'transient-heat')
  exec_3d = os.path.join(this_dir, '..', '..', build_folder, 'examples', 'TSHT', 'transient-heat')
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

base_cfg = {'dim': 3, "BDF2": "false", 'forcing': 'true', }

cases = [
    # check serial/parallel for linear/quadratic, and mfree/matrix
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "false", 'ntasks': 1},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "false", 'ntasks': 1},
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "true", 'ntasks': 1},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "true", 'ntasks': 1},  # quite slow
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "false", 'ntasks': 4},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "false", 'ntasks': 4},
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "true", 'ntasks': 4},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'mfree': "true", 'ntasks': 4},
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'forcing': 'false', 'mfree': "false", 'ntasks': 4},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'forcing': 'false', 'mfree': "false", 'ntasks': 4},
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.1, 'forcing': 'false', 'mfree': "true", 'ntasks': 4},
    {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.1, 'forcing': 'false', 'mfree': "true", 'ntasks': 4},

    # user can check the slope of linear basis functions (not included in the regression test)
    {'order': "linear", 'refine_lvl': 2, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 2},
    {'order': "linear", 'refine_lvl': 3, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 4},
    {'order': "linear", 'refine_lvl': 4, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 4},
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 8},

    # BDF2 & BE timestepper
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.2, 'BDF2': "false", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.1, 'BDF2': "false", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.05, 'BDF2': "false", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},

    {'order': "linear", 'refine_lvl': 5, 'dt': 0.2, 'BDF2': "true", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.1, 'BDF2': "true", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},
    {'order': "linear", 'refine_lvl': 5, 'dt': 0.05, 'BDF2': "true", 'forcing': 'true', 'mfree': "false", 'ntasks': 8},


    # user can check the slope of quadratic basis functions (not included in the regression test)
    # {'order': "quadratic", 'refine_lvl': 2, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 2},
    # {'order': "quadratic", 'refine_lvl': 3, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 4},
    # {'order': "quadratic", 'refine_lvl': 4, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 4},
    # {'order': "quadratic", 'refine_lvl': 5, 'dt': 0.05, 'BDF2': "true", 'mfree': "false", 'ntasks': 8},
]

from regression import RegressionTester, RegexDiffMetric, VecDiffMetric

reg = RegressionTester()

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
reg.add_metric(RegexDiffMetric('Final L2 error', 'output.txt', re.compile(r'L2Error : (.*)'), 5e-8))
reg.add_file_generator("config.txt", """
basisFunction = "{order}"
mfree = {mfree}

dt = {dt}
totalT = 1.0
BDF2 = {BDF2}
forcing = {forcing}
OutputStartTime = 0  # do not output anything other than initial condition
OutputInterval = 1

background_mesh = {{
  refine_lvl = {refine_lvl}
  min = [0, 0, 0]
  max = [1, 1, 1]
}}

dump_vec = true

#################### solver setting ####################
solver_options_ht = {{
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "asm"
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

reg.add_cases(base_cfg, cases)
reg.set_exclude_patterns(["*.vtu", "*.pvtu", "*.info"])

reg.add_metric(VecDiffMetric('solution_vec.vec', 5e-8))
reg.add_metric(RegexDiffMetric('Final L2 error', 'output.txt', re.compile(r'L2Error : (.*)'), 5e-8))
reg.add_file_generator("config.txt", """
basisFunction = "{order}"
mfree = {mfree}

dt = {dt}
totalT = 1.0
BDF2 = {BDF2}
forcing = {forcing}
OutputStartTime = 0  # do not output anything other than initial condition
OutputInterval = 1

background_mesh = {{
  refine_lvl = {refine_lvl}
  min = [0, 0, 0]
  max = [1, 1, 1]
}}

dump_vec = true

#################### solver setting ####################
solver_options_ht = {{
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "asm"
  ksp_atol = 1e-7
  ksp_rtol = 1e-10
  ksp_converged_reason = ""
  ksp_monitor = ""
}}
""")
reg.set_run_cmd(executable_2d)
reg.run_main()
