#!/usr/bin/env python
import os
import sys
import shutil
import subprocess

sys.path.append(os.environ['TALYFEM_DIR'] + "/python_scripts")
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS

config_file = "config.txt"
EXECUTABLE_PATH = os.environ['EXECUTABLE_PATH']

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def write_config(dd_val):
  true_false = ["False", "True"]
  config = open(config_file, "w")
  config.write("nsd = 2\n")
  config.write("Lx = 8\n")
  config.write("Nelemx = 4\n")
  config.write("Ly = 8\n")
  config.write("Nelemy = 4\n")
  config.write("ifBoxGrid = True\n")
  config.write("typeOfIC = 0\n")
  config.write("orderOfBF = 1\n")
  config.write("ifDD = " + true_false[dd_val] +"\n")


def do_test():
  res = 0

  # serial test
  for dd_val in (0, 1):
    write_config(dd_val)
    comm = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', '1', EXECUTABLE_PATH]
    null = open(os.devnull, "w")
    print comm
    res += subprocess.call(comm, stdout=null, stderr=null)
    null.close()

  # parallel test
  for dd_val in (0, 1):
    write_config(dd_val)
    comm = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', '3', EXECUTABLE_PATH]
    null = open(os.devnull, "w")
    print comm
    res += subprocess.call(comm, stdout=null, stderr=null)
    null.close()

  if res == 0:
    print "...test " + PASS_COLOR + "passed" + END_COLOR
  else:
    print "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
  return res


if __name__ == '__main__':
  rval = do_test()
  sys.exit(rval)
