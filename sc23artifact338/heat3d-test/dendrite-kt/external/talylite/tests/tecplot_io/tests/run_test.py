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
PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def do_test(program, filename_out, filename_err):
  shutil.copyfile(TEST_DATA_DIR + '/data.plt', 'data.plt')

  sout = open(filename_out, "w")
  serr = open(filename_err, "w")
  subprocess.call(program, stdout=sout, stderr=serr)
  sout.close()
  serr.close()

  filter_file(filename_out)
  filter_file(filename_err)

  comm = [
    ["no warnings", ["diff", filename_out, TEST_DATA_DIR + "/expected_" + filename_out]], 
    ["no errors", ["diff", filename_err,  TEST_DATA_DIR + "/expected_" + filename_err]],
    ["output is the same as input", ["diff", TEST_DATA_DIR + "/data.plt", "data_copy.plt"]],
    ["manual data output is correct", ["diff", "manual_data.plt", TEST_DATA_DIR + "/manual_data_expected.plt"]],
    ["multi-zone manual data output is correct", ["diff", "manual_data_multizone.plt", TEST_DATA_DIR + "/expected_manual_data_multizone.plt"]]
    ]

  null = open(os.devnull, "w")
  error = 0
  for value in comm:
    res = subprocess.call(value[1], stdout=null)
    if res == 0:
      print value[0] + "...test " + PASS_COLOR + "passed" + END_COLOR
    else:
      print value[0] + "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
      error = 1

  null.close()
  return error


if __name__ == '__main__':
  import sys

  null = open(os.devnull, "w")
  subprocess.call(["rm", "-f", "data_copy.plt"], stdout=null)
  subprocess.call(["rm", "-f", "manual_data_multizone.plt"], stdout=null)
  null.close()

  print "testing basic tecplotIO serial program",

  serial_program = [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ['-n', '1', EXECUTABLE_PATH]
  rval = do_test(serial_program, "out.txt", "err.txt")

  sys.exit(rval)
