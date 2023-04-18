#!/usr/bin/env python
import os
import shutil
import subprocess
import sys

sys.path.append(os.environ['TALYFEM_DIR'] + "/python_scripts")
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS
from test_common import filter_file

EXECUTABLE_PATH = os.environ['EXECUTABLE_PATH']
TEST_DATA_DIR = os.environ['TEST_DATA_DIR']
program = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ["-n", "1", EXECUTABLE_PATH]
config_file = "config.txt"
config_file_base = TEST_DATA_DIR + "/config_base.txt"

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def set_config_file(nsd, nDoF, perBounds, perVars):
  shutil.copyfile(config_file_base, config_file)
  with open(config_file, 'a') as f:
    f.write("nsd = " + str(nsd) + "\n")
    f.write("nDoF = " + str(nDoF) + "\n")
    f.write("nPerBounds = " + str(len(perBounds)) + "\n")
    f.write("nPerVars = " + str(len(perVars)) + "\n")
    for i, bound in enumerate(perBounds):
      name = "perBounds" + str(i+1)
      f.write(name + " = " + str(bound) + "\n")
    for i, var in enumerate(perVars):
      name = "perVars" + str(i+1)
      f.write(name + " = " + str(var) + "\n")


def do_run(output_file, detail_line):
  use_valgrind = False
  program_prefix = []
  if use_valgrind:
    program_prefix = ["valgrind", "--tool=memcheck", "--leak-check=full", 
                      "-v", "--show-reachable=yes"]
  with open(output_file, 'a') as f:
    f.write(detail_line + ":\n")
    f.flush()
    if use_valgrind:
      subprocess.call(program_prefix + program)
    else:
      subprocess.call(program_prefix + program, stdout=f)


def do_test(filename, data_values):
  error = 0
  with open(os.devnull, "w") as null:
    print "testing periodic fill", filename, 
    with open(filename, 'w') as f:
      f.write("")
    for (nsd, ndof, perbounds, pervars) in data_values:
      set_config_file(nsd, ndof, perbounds, pervars)
      details = ("nsd = " + str(nsd) + " ndof = " + str(ndof) +
                 " periodicBounds = " + str(perbounds) +
                 " periodicVariables = " +  str(pervars))
      do_run(filename, details)

    filter_file(filename)
    comm = ["diff", filename,  TEST_DATA_DIR + "/expected_" + filename]
    res = subprocess.call(comm, stdout=null)
    if res == 0:
      print "...test " + PASS_COLOR + "passed" + END_COLOR
    else:
      print "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
      error = 1
  return error


if __name__ == '__main__':
  outputfile = "output.txt"
  input_data = [
   # nsd, nDoF, perBounds, perVars
   [1, 1, [2], [0]],
   [1, 4, [2], [2]],
   [1, 4, [2], [0, 1, 2]],
   [2, 2, [4], [1, 0]],
   [2, 2, [2, 4], [1, 0]],
   [3, 3, [2, 4, 6], [0, 2, 1]],
   [3, 2, [6], [1, 0]],
   [3, 1, [2, 4], [0]],
   [3, 5, [2, 6], [2, 0, 4]],
   [3, 5, [2, 4, 6], [2]],
  ] 
  rval = do_test(outputfile, input_data)

  import sys
  sys.exit(rval)
