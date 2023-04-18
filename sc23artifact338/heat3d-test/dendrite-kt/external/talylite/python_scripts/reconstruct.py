#!/usr/bin/python

import argparse
import os
import stat
import sys
import subprocess
import fnmatch
import base64

# Load libconf from local directory
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/external")
import libconf

def in_unix_list(val, lst):  
  for v in lst:
    if fnmatch.fnmatch(val, v):
      return True
  return False

def rebuild_code(build_info, output_dir, source_repo):
  # gather config variables early to fail fast
  start_hash = build_info['git']['last_commit_hash']
  outstanding_changes = (build_info['git']['uncommitted_changes'] != '0')
  patch = build_info['git']['patch']

  old_cwd = os.getcwd()
  try:
    # make sure we're not going to clobber an existing directory
    if os.path.exists(output_dir):
      raise Exception('Output directory "' + output_dir + '" already exists')

    # clone git repo
    try:
      os.makedirs(output_dir)
      print source_repo, output_dir
      subprocess.check_call(['git', 'clone', source_repo, output_dir])
    except subprocess.CalledProcessError:
      raise Exception('Unable to clone Git directory "' + str(source_repo) + '" into "' + str(output_dir) + '"')

    os.chdir(output_dir)

    # check out last commit
    try:
      subprocess.check_call(['git', 'checkout', start_hash])
    except subprocess.CalledProcessError:
      raise Exception('Unable to checkout hash "' + start_hash + '" - repository "' + str(source_repo) + "' likely does not contain the required commit.")

    # apply uncommitted changes
    if outstanding_changes:
      with open('patch.diff', 'w') as f:
        f.write(patch)
      subprocess.check_call(['git', 'apply', 'patch.diff'])
      #os.remove('patch.diff')

  finally:
    os.chdir(old_cwd)

def rebuild_runtime(runtime_info, output_dir, include_run_env_vars):
  # gather variables early to fail fast
  cmd_line_args = runtime_info['cmd_line_args']
  n_procs = runtime_info['n_procs']
  run_env = dict(runtime_info['environment']) # convert list of tuples to dict
  input_files = runtime_info.get('input_files', [])

  old_cwd = os.getcwd()
  try:
    if os.path.exists(output_dir):
      raise Exception('Output directory "' + output_dir + '" already exists')

    # make run directory
    os.makedirs(output_dir)
    os.chdir(output_dir)

    # write run environment
    with open('environment.sh', 'w') as f:
      for key, value in sorted(run_env.iteritems()):
        if not in_unix_list(key, include_run_env_vars):
          f.write('#')
        f.write(key + '=' + value + '\n')

    # write run script
    with open('run.sh', 'w') as f:
      run_prefix = 'mpirun -n ' + n_procs
      executable = '../source/project/' + cmd_line_args[0].split('/')[-1]
      f.write('# This is an automatically generated file.\n\n')
      f.write('source ./environment.sh\n\n')
      f.write(run_prefix + ' ' + executable + ' ' + (' '.join(cmd_line_args[1:])))

    # write input files
    for pair in input_files:
      name = pair[0]
      with open(name, 'w') as f:
        f.write(base64.b64decode(pair[1]))

    # mark run script as user executable
    st = os.stat('run.sh')
    os.chmod('run.sh', st.st_mode | stat.S_IEXEC)

  finally:
    os.chdir(old_cwd)

CFG_BUILDINFO_TALY = 'BuildInfo_taly'
#CFG_BUILDINFO_PROJECT = 'BuildInfo_taly_project'
CFG_RUNTIME = 'runtime'

parser = argparse.ArgumentParser(description='Reconstruct a TalyFEM project from results.')

parser.add_argument('--taly-git', default='https://bitbucket.org/baskargroup/taly_fem', type=str,
                    help='Git repository path to consider pulling from when rebuilding source code. May be local paths.')
#parser.add_argument('--project-git', type=str, default=None,
#                    help='Git repository path to consider pulling from when rebuilding source code. May be local paths.')
parser.add_argument('--run-out', type=str, default="run", help='Output directory for run working directory.')
parser.add_argument('--run-env-vars', type=str, default=['PETSC_*', 'LD_LIBRARY_PATH', 'PATH'],
                    help='List of environment variables to include in run environment. Wildcards are supported.')
parser.add_argument('repro_file', default='repro.cfg', type=argparse.FileType('r'), help='Path to TalyFEM repro.cfg.')

args = parser.parse_args()

config = libconf.load(args.repro_file)

rebuild_code(config[CFG_BUILDINFO_TALY], 'source/taly', args.taly_git)
rebuild_runtime(config[CFG_RUNTIME], args.run_out, args.run_env_vars)
