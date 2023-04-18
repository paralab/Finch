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
program = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ["-n", "2", EXECUTABLE_PATH]
config_file = "config.txt"
config_file_base = TEST_DATA_DIR + "/config_base.txt"

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def set_config_file(nsd, expected_measure, if_box_grid, input_grid_filename):
  shutil.copyfile(config_file_base, config_file)
  with open(config_file, 'a') as f:
    f.write("nsd = " + str(nsd) + "\n")
    f.write("ifBoxGrid = " + str(if_box_grid) + "\n")
    f.write("expected_measure = " + str(expected_measure) + "\n")
    if not if_box_grid:
      f.write('inputFilenameGrid = "' + str(input_grid_filename) + '"\n')


def do_run(output_file, detail_line):
  program_prefix = []
  with open(output_file, 'a') as f:
    f.write(detail_line + ":\n")
    f.flush()
    res = subprocess.call(program_prefix + program, stdout=f, stderr=f)
    f.write("\n")
    return res


def do_test(data_values):
  error = 0
  with open(os.devnull, "w") as null:
    print "Running Integrator Test:"
    with open("output.txt", 'w') as f:
      f.write("")
    for (testname, nsd, expected_measure, if_box_grid,
         input_grid_filename) in data_values:
      print "testing ", testname,
      set_config_file(nsd, expected_measure, if_box_grid, input_grid_filename)
      details = ("nsd = " + str(nsd) + " ifBoxGrid = " + str(if_box_grid) +
                 " inputFilenameGrid = " + str(input_grid_filename))
      res = do_run("output.txt", details)

      if res == 0:
        print "...test " + PASS_COLOR + "passed" + END_COLOR
      else:
        print "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
        error = 1
  return error


if __name__ == '__main__':
  input_data = [
   # (testname, nsd, expected_measure, ifBoxGrid, filename)
   ("1D box grid", 1, 1.0, True, ""),
   ("2D box grid", 2, 1.0, True, ""),
   ("3D box grid", 3, 1.0, True, ""),
  ]
  rval = do_test(input_data)

  import sys
  sys.exit(rval)
