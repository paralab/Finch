#!/usr/bin/env python
import os
import sys
import shutil
import subprocess

sys.path.append(os.environ['TALYFEM_DIR'] + "/python_scripts")
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS

EXECUTABLE_PATH = os.environ['EXECUTABLE_PATH']
TEST_DATA_DIR = os.environ['TEST_DATA_DIR']
config_file = "config.txt"
config_file_base = TEST_DATA_DIR + "/config_base.txt"

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def set_config_file(boundaries, variables):
  shutil.copyfile(config_file_base, config_file)
  with open(config_file, 'a') as f:
    dof = len(variables)
    f.write("dof = " + str(dof) + "\n")
    for i, val in enumerate(variables):
      if val:
        f.write("isPerVar" + str(i+1) + " = True\n")
      else:
        f.write("isPerVar" + str(i+1) + " = False\n")
    for key, val in zip(['X', 'Y', 'Z'], boundaries):
      f.write("boundaryType" + key + " = \"" + val + "\"\n")


def do_test(n_procs, output_filename, test_cases):
  res = 0
  # check for memory issues, this needs to be done manually
  use_valgrind = False
  program_prefix = []
  if use_valgrind:
    program_prefix = ["valgrind", "--tool=memcheck", "--leak-check=full", 
                      "-v", "--show-reachable=yes"]

  program = program_prefix + [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', str(n_procs), EXECUTABLE_PATH]
  print program

  f = open(output_filename, 'w')
  all_passed = True
  with open(os.devnull, "w") as null:
    for boundaries, variables in test_cases:
      set_config_file(boundaries, variables)
      f.write("\n==========================================================\n")
      f.write("==========================================================\n")
      for key, val in zip(['X', 'Y', 'Z'], boundaries):
        f.write(" boundaryType" + key + " = " + val)
      f.write("\n")
      for i, val in enumerate(variables):
        f.write(" isPerVar" + str(i+1) +  " = " + str(val))
      f.write("\n")
      f.flush()
      if use_valgrind:
        if subprocess.call(program) != 0:
          all_passed = False
      else:
        if subprocess.call(program, stdout=null, stderr=f) != 0:
          all_passed = False
      f.flush()

    if all_passed:
      print "...test " + PASS_COLOR + "passed" + END_COLOR
    else:
      print "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
  f.close()
  if not all_passed:
    res = 1
  return res


if __name__ == '__main__':
  rval = 0
  for n_procs in (1, 2, 3):
    outfile = "output_np" + str(n_procs) + ".txt"
    
    test_cases = []
    # no periodic boundaries
    boundary_values = ['direchlet', 'direchlet']
    per_var = [False]
    test_cases.append([boundary_values, per_var])
    per_var = [False, False]
    test_cases.append([boundary_values, per_var])
    # periodic in x
    boundary_values = ['periodic', 'direchlet']
    per_var = [True]
    test_cases.append([boundary_values, per_var])
    for perx in (True, False):
      for pery in (True, False):
        if (not perx) and (not perx): 
          continue
        per_var = [perx, pery]
        test_cases.append([boundary_values, per_var])
    # periodic in y
    boundary_values = ['direchlet', 'periodic']
    per_var = [True]
    test_cases.append([boundary_values, per_var])
    for perx in (True, False):
      for pery in (True, False):
        if (not perx) and (not perx): 
          continue
        per_var = [perx, pery]
        test_cases.append([boundary_values, per_var])
    # periodic in x, y
    boundary_values = ['periodic', 'periodic']
    per_var = [True]
    test_cases.append([boundary_values, per_var])
    for perx in (True, False):
      for pery in (True, False):
        if (not perx) and (not perx): 
          continue
        per_var = [perx, pery]
        test_cases.append([boundary_values, per_var])
    print "testing with np =", n_procs,
    rval += do_test(n_procs, outfile, test_cases)
  sys.exit(rval)
