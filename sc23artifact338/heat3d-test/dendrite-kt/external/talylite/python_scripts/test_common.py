#!/usr/bin/python

import os
import sys
import shutil
import subprocess
import stat
import time

# libconf requires Python 2.7 or above, or we get syntax errors
if sys.version_info[1] < 7:
  print("Python 2.7 or higher required. Try 'module load python/2.7'?")
  sys.exit(1)

# Load libconf from local directory
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/external")
import libconf

# Runs a series of tests.
# A "test report" will be automatically printed after all tests complete
# (or are interrupted with Ctrl-C).

# Usage:

# from test_common import run_tests
# run_tests(test_list, program_path, program_args)
# OR
# run_tests(test_list, program_path, program_args, result_processing_func)

# Where test_list is a list of dictionaries representing test config.txts, 
# program_path is the path to the executable to launch, and
# program_args is an array of command-line arguments.

# IMPORTANT: YOU SHOULD USE ABSOLUTE PATHS IN BOTH PROGRAM_PATH AND ARGUMENTS!
# Since the test will be run from results/[test_id]/, relative paths 
# usually need to be prepended with "../../".
# Use os.getcwd() + "/my_program" if you want to turn a relative path to 
# into an absolute one.

# Optionally, you can also supply a result_processing_func if your program
# does not follow the Tutorials test output format.

# If you want to write your own, the processing function looks like:
# process(test, test_result_directory)
# Where test is the usual dictionary format and test_result_directory is the
# directory the test was run in. Program stderr and stdout is redirected to
# test_result_directory/OUT.LOG, so you probably want to parse that.
# You should return a tuple in the form of (boolean_passed, msgs_list), where
# boolean_result is True if the test passed or False if the test failed,
# and msgs_list is a list of messages (only displayed if the test failed).
# If you want to print out result information as tests complete, just print
# it directly in this function.

# The test format is a simple dictionary, like so:

# non_dd = {
#   "__DESCR": "An example test.",
#   "nsd": 2,
#   "ifBoxGrid": 1,
#   "orderOfBF": 1,
#   "ifDD": 0,
#   "Lx": 1,
#   "Nelemx": 20,
#   "Ly": 1,
#   "Nelemy": 20,
#   "typeOfIC": 1,
#   "dt": 0.001,
#   "nOfTS": 10,
#   "outputExtension": ".plt"
# }

# Any keys that begin with "__" will be IGNORED by the config file writer and
# parameter expansion.

# Tests are run in "./results/test_[test_idx]".
# Program output is redirected to an OUT.LOG file in the same directory.

# A test description can be set with the "__DESCR" key/value pair.
# It will be written as a comment in the top of the test's config.txt.

# A test can be run with extra arguments with a "__PROGRAM_ARGS" key/value
# pair. The value type should be list. The arguments will be appended to
# program_args.

# If a test requires other files in its directory to work, you can use
# "__COPY_FILES": ["file1.plt", "file2.plt"]
# to have them automatically copied to the test directory.

# If a test should be *skipped*, you can add "__SKIP" = "reason" to the
# test. You can use this instead of removing the test from the list to
# provide a warning that the test was skipped and keep consistent test IDs.

# You can also use "__" variables as part of command line arguments, and they
# will be automatically replaced if set for that test. 
# For example:
#   test = {
#    "__NUM_PROCESSES": 2
#   }
#   run_tests([test], "mpirun", ["-np", "__NUM_PROCESSES", "CODE/HT"])
# executes the program as:
#   mpirun -np 2 CODE/HT

PASS_COLOR = '\033[92m'
FAIL_COLOR = '\033[91m'
END_COLOR = '\033[0m'

def write_config(path, test):
  with open(path, 'w') as f:
    f.write("# THIS IS AN AUTOMATICALLY GENERATED FILE.\n")

    if test.has_key("__DESCR"):
      f.write("# [DESCR]: " + test["__DESCR"] + "\n")

    # strip entries that start with _
    test_clean = {}
    for key, val in test.iteritems():
      if not key.startswith("__"):
        test_clean[key] = val

    # dump to file
    libconf.dump(test_clean, f)

def run_test(test, test_id, program_path, program_args):
  test_dir = "results/test_" + str(test_id)

  # clean the directory
  if os.path.exists(test_dir):
    shutil.rmtree(test_dir)

  os.makedirs(test_dir)

  # copy files
  if test.has_key("__COPY_FILES"):
    files = test["__COPY_FILES"]
    if type(files) is str:
      files = [test["__COPY_FILES"]]

    for f in files:
      # if f doesn't exist as-is, try the TEST_DATA_DIR environment variable
      if not os.path.exists(f):
        if "TEST_DATA_DIR" in os.environ:
          data_dir_f = os.environ["TEST_DATA_DIR"] + "/" + f
          if os.path.exists(data_dir_f):
            f = data_dir_f

      if os.path.isdir(f):
        shutil.copytree(f, test_dir)
      else:
        shutil.copy(f, test_dir)

  write_config(test_dir + "/config.txt", test)
  test_output_path = test_dir + "/OUT.LOG"

  # replace args in the form of "__VAR" with test["__VAR"]
  call = [program_path]
  for arg in program_args:
    if arg.startswith("__") and test.has_key(arg):
      call.append(str(test[arg]))
    else:
      call.append(os.path.expandvars(arg))

  # write a run script, in case you want to run the test yourself
  run_script_path = test_dir + "/run_test.sh"
  with open(run_script_path, 'w') as run_script:
    run_script.write("#!/bin/bash\n\n")
    run_script.write("# THIS IS AN AUTOMATICALLY GENERATED FILE.\n")
    run_script.write(" ".join(call))
  # make it executable
  os.chmod(run_script_path, os.stat(run_script_path).st_mode | stat.S_IXUSR)

  # run the test
  start_time = time.time()
  with open(test_output_path, 'w') as output_file:
    subprocess.call(call, stdout=output_file, stderr=output_file, cwd=test_dir)
  end_time = time.time()
  print "  (done, in " + str(round(end_time - start_time, 2)) + "s)"

  return test_dir

def process_default(test, test_dir):
  result = False
  error = False
  msgs = []
  with open(test_dir + "/OUT.LOG") as f:
    for line in f:
      if "PETSC ERROR" in line and not error:
        error = True
        msgs.append(FAIL_COLOR + "PETSC ERROR(S) DETECTED" + END_COLOR)
        print "  " + FAIL_COLOR + "PETSC ERROR(S) DETECTED" + END_COLOR

      if "TEST" in line or "ERR" in line or "WARNING" in line or "error" in line:
        if "PETSC ERROR" not in line:
          msgs.append(line.strip())

      if PASS_COLOR in line and not error:
        print "  " + line.strip()
        result = True
      elif FAIL_COLOR in line:
        print "   " + line.strip()

  return (result and not error, msgs)

def print_report(results):
  total = len(results)
  passed = sum(result[0] == True for result in results)
  failed = total - passed

  print ""
  print "TEST REPORT"
  print "==========="
  if passed == total:
    print PASS_COLOR + str(passed) + "/" + str(total) + " tests passed!" + END_COLOR
  else:
    print str(passed) + "/" + str(total) + " tests passed!"

  for i in range(0, total):
    result = results[i]
    test_id = result[2]
    if result[0] == False:
      print FAIL_COLOR + "Test " + str(test_id) + " failed!" + END_COLOR
      for msg in result[1]:
        print "   " + msg

def run_tests(tests, program_path, program_args = [], process_result_func = process_default):
  # results is a list of tuples in the form of
  # (test_passed, [msg1, msg2, ...], test_id)
  results = []

  test_ids_to_run = []
  for arg in sys.argv:
    if arg.isdigit():
      test_ids_to_run.append(int(arg))

  # if the user didn't supply a list of tests to run, run everything
  if len(test_ids_to_run) == 0:
    test_ids_to_run = range(0, len(tests))

  try:
    for i in test_ids_to_run:
      test_args = list(program_args)
      if tests[i].has_key("__SKIP"):
        print "SKIPPING TEST " + str(i) + ": " + tests[i]["__SKIP"]
        continue

      # plain test, run it directly
      if tests[i].has_key("__PROGRAM_ARGS"):
        test_args += tests[i]["__PROGRAM_ARGS"]

      test_id = str(i)
      perc_done = int(round(test_ids_to_run.index(i) / float(len(test_ids_to_run)) * 100))
      print "Running test " + str(i) + "... (" + str(perc_done) + "%)"

      # actually run the test
      test_output_path = run_test(tests[i], test_id, program_path, test_args)
      if process_result_func:
        result = process_result_func(tests[i], test_output_path)
        results.append(result + (test_id,))

  except KeyboardInterrupt:
    print "Testing cancelled."

  print_report(results)
  return results

DEFAULT_PROGRAM_PATH = os.getenv("TALYFEM_RUN_COMMAND")
if DEFAULT_PROGRAM_PATH == None:
  DEFAULT_PROGRAM_PATH = "mpirun"

DEFAULT_PROGRAM_ARGS = []
if DEFAULT_PROGRAM_PATH == "aprun":
  DEFAULT_PROGRAM_ARGS += ['-q']
if DEFAULT_PROGRAM_PATH == "ibrun":
  DEFAULT_PROGRAM_ARGS += ['-o', '0']

def have_mumps():
  path = os.getenv("PETSC_DIR") + os.sep + os.getenv("PETSC_ARCH", "") + os.sep + "lib" + os.sep + "petsc" + os.sep + "conf" + os.sep + "petscvariables"
  found_mumps = False
  with open(path, "r") as f:
    for line in f:
      if line.startswith("TEST_RUNS ="):
        found_mumps = ("MUMPS" in line)
  return found_mumps

def filter_file(path):
  if DEFAULT_PROGRAM_PATH == "ibrun":
    lines = [line for line in open(path) if (not line.startswith("TACC:") and not line == " \n")]
    open(path, 'w').writelines(lines)

def all_tests_passed(results):
  for result in results:
    if result[0] == False:
      return False
  return True

