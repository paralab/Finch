#!/usr/bin/env python
import os
import shutil
import subprocess
import sys

sys.path.append(os.getcwd() + "/../../../PythonScripts")
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS
from test_common import filter_file

config_file = "config.txt"
config_file_base = "config_base.txt"

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def set_config_file(values, filename, domain_decomp):
  shutil.copyfile(config_file_base, config_file)
  with open(config_file, 'a') as f:
    dims = len(values)
    f.write("nsd = " + str(dims) + "\n")
    f.write("outputFile = \"" + filename + "\"\n")
    if domain_decomp:
      f.write("ifDD = 1\n")
    else:
      f.write("ifDD = 0\n")
    for key, val in zip(['X', 'Y', 'Z'], values):
      f.write("boundaryType" + key + " = \"" + val + "\"\n")


def do_test(program, n_dims, filename_base, boundary_values, domain_decomp):
  err = 0
  res = 0
  use_valgrind = False
  program_prefix = []
  if use_valgrind:
    program_prefix = ["valgrind", "--tool=memcheck", "--leak-check=full", 
                      "-v", "--show-reachable=yes", "--track-origins=yes"]

  diff_dir_prefix = ""
  if os.getenv("ON_BW"):
    diff_dir_prefix = "bw/"

  with open(os.devnull, "w") as null:
    for values in boundary_values:
      print "  testing", "_".join(values), "..."
      filename = filename_base + str(n_dims) + "_" + "_".join(values) + ".plt"
      set_config_file(values, filename, domain_decomp)
      if use_valgrind:
        subprocess.call(program_prefix + program)
      else:
        subprocess.call(program_prefix + program, stdout=null)

      filter_file(filename)
      comm = ["diff", '-w', filename,  diff_dir_prefix + "expected_" + filename]
      res = subprocess.call(comm) #, stdout=null)

      if res == 0:
        print "  ...test " + PASS_COLOR + "passed" + END_COLOR
      else:
        print "  ...test " + FAIL_COLOR + "***FAILED***" + END_COLOR
        print "     (" + " ".join(comm) + ")"
        err = 1

  return err


if __name__ == '__main__':
  boundary_values = [0,0,0,0]
  boundary_values[1] = [
   ['periodic'],
   ['direchlet'],
  ]
  boundary_values[2] = [
   ['periodic', 'periodic'],
   ['direchlet', 'periodic'],
   ['periodic', 'direchlet'],
   ['direchlet', 'direchlet'],
  ]
  boundary_values[3] = [
   ['periodic', 'periodic', 'periodic'],
   ['direchlet', 'periodic', 'periodic'],
   ['periodic', 'direchlet', 'periodic'],
   ['direchlet', 'direchlet', 'periodic'],
   ['periodic', 'periodic', 'direchlet'],
   ['direchlet', 'periodic', 'direchlet'],
   ['periodic', 'direchlet', 'direchlet'],
   ['direchlet', 'direchlet', 'direchlet'],
  ]

  rval = 0
  executable = '../pbc_node_test'
  program = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', '1', executable]
  for n_dims in (1,2,3):
    print "testing", n_dims, "dimension serial... "
    rval += do_test(program, n_dims, "out", boundary_values[n_dims], False)

  for n_dims in (2, 3):
    print "testing", n_dims, "dimension serial, with domain decomp... "
    rval += do_test(program, n_dims, "dd1_out", boundary_values[n_dims], True)
  
  program = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', '3', executable]
  for n_dims in (2, 3):
    print "testing", n_dims, "dimension parallel, no domain decomp... "
    rval += do_test(program, n_dims, "out", boundary_values[n_dims], False)

  for n_dims in (2, 3):
    print "testing", n_dims, "dimension parallel, with domain decomp... "
    rval += do_test(program, n_dims, "dd3_out", boundary_values[n_dims], True)

  import sys
  sys.exit(rval)
