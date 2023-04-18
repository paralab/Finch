#!/usr/bin/env python
import os
import sys
import shutil
import subprocess

sys.path.append(os.environ['TALYFEM_DIR'] + "/python_scripts")
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS
from test_common import filter_file

EXECUTABLE_PATH = os.environ['EXECUTABLE_PATH']
TEST_DATA_DIR = os.environ['TEST_DATA_DIR']
config_file = "config.txt"
config_file_base = TEST_DATA_DIR + "/config_base.txt"

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def set_config_file(values):
  shutil.copyfile(config_file_base, config_file)
  with open(config_file, 'a') as f:
    for key, val in values.items():
      if val:
        f.write(key + " = True\n")
      else:
        f.write(key + " = False\n")


def do_test(program, filename_out, filename_err, data):
  sout = open(filename_out, "w")
  serr = open(filename_err, "w")
  for description, values in data:
    set_config_file(values)
    sout.write("*****" + description + "\n");
    sout.flush();
    serr.write("*****" + description + "\n");
    serr.flush();
    subprocess.call(program, stdout=sout, stderr=serr)
  sout.close()
  serr.close()

  # remove "TACC:" echoes so diff doesnt fail
  filter_file(filename_out)
  filter_file(filename_err)

  comm1 = ["diff", filename_out,  TEST_DATA_DIR + "/expected_" + filename_out]
  comm2 = ["diff", filename_err,  TEST_DATA_DIR + "/expected_" + filename_err]

  null = open(os.devnull, "w")
  res = subprocess.call(comm1, stdout=null)
  res += subprocess.call(comm2, stdout=null)
  null.close()

  if res == 0:
    print "...test " + PASS_COLOR + "passed" + END_COLOR
  else:
    print "...test " + FAIL_COLOR + "***FAIILED***" + END_COLOR
  return res


if __name__ == '__main__':
  import sys
  
  ERR_KEY = "ifPrintError"
  INFO_KEY = "ifPrintInfo"
  LOG_KEY = "ifPrintLog"
  STAT_KEY = "ifPrintStat"
  WARN_KEY = "ifPrintWarn"

  data = []
  
  desc = "Using default values for all (should all print)"
  data.append([desc, {}])

  desc = "Turning off ERROR (should have no effect and be same as above)"
  values = {ERR_KEY: False}
  data.append([desc, values])

  desc = "Turning on everything (should be same as above, except input echo)"
  values = {ERR_KEY: True, INFO_KEY: True, LOG_KEY: True, 
            STAT_KEY: True, WARN_KEY: True}
  data.append([desc, values])

  desc = "Turning off INFO"
  values = {INFO_KEY: False}
  data.append([desc, values])

  desc = "Turning off LOG"
  values = {LOG_KEY: False}
  data.append([desc, values])

  desc = "Turning off STAT"
  values = {STAT_KEY: False}
  data.append([desc, values])

  desc = "Turning off WARN"
  values = {WARN_KEY: False}
  data.append([desc, values])

  desc = "Turning off everything (only ERROR and (cout)/(cerr) should print)"
  values = {ERR_KEY: False, INFO_KEY: False, LOG_KEY: False, 
            STAT_KEY: False, WARN_KEY: False}
  data.append([desc, values])

  # check for memory issues, this needs to be done manually
  use_valgrind = False
  program_prefix = []
  if use_valgrind:
    program_prefix = ["valgrind", "--tool=memcheck", "--leak-check=full", 
                      "--show-reachable=yes"]

  rval = 0

  print "testing serial program",
  serial_program = program_prefix + [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ["-n", "1", EXECUTABLE_PATH]
  rval += do_test(serial_program, "out.txt", "err.txt", data)

  print "testing parallel program",
  parallel_program = program_prefix + [DEFAULT_PROGRAM_PATH] + DEFAULT_PROGRAM_ARGS + ["-n", "2", EXECUTABLE_PATH]
  rval += do_test(parallel_program, "out.txt", "err.txt", data)
  
  sys.exit(rval)
