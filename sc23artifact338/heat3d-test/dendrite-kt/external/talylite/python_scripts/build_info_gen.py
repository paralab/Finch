#!/usr/bin/env python

import socket
import os
import sys
import traceback
import calendar
import subprocess
import getpass
from datetime import datetime

# Load argparse from local directory if it isn't installed.
# We run on some Python 2.6 systems that either:
#  (1) don't have the argparse module,
#  (2) have a buggy argparse module that always creates build_info.h
#      even when --output is specified
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/external")
import argparse

# Gonna need a shower after this one
# http://stackoverflow.com/a/34055481
if not hasattr(subprocess, 'check_output'):
    # python 2.6 doesn't include check_output
    # monkey patch it in!
    import subprocess
    STDOUT = subprocess.STDOUT

    def check_output(*popenargs, **kwargs):
        if 'stdout' in kwargs:  # pragma: no cover
            raise ValueError('stdout argument not allowed, '
                             'it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE,
                                   *popenargs, **kwargs)
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd,
                                                output=output)
        return output
    subprocess.check_output = check_output

    # overwrite CalledProcessError due to `output`
    # keyword not being available (in 2.6)
    class CalledProcessError(Exception):

        def __init__(self, returncode, cmd, output=None):
            self.returncode = returncode
            self.cmd = cmd
            self.output = output

        def __str__(self):
            return "Command '%s' returned non-zero exit status %d" % (
                self.cmd, self.returncode)
    subprocess.CalledProcessError = CalledProcessError




DEVNULL = open(os.devnull, 'w')

# register("gen/timestamp_str", str) 
#  =>
# categories["gen"]["timestamp_str"] == (str, timestamp_str)
categories = {}

class rawstr:
  pass
class uint64_t:
  pass

def register(path, type):
  def wrapper(f):
    path_arr = path.split('/')
    category = path_arr[0]
    name = path_arr[1]

    category_dict = categories.get(category, {})
    category_dict[name] = (type, f)
    categories[category] = category_dict
    return f
  return wrapper


@register("gen/timestamp", uint64_t)
def timestamp():
  return calendar.timegm(datetime.utcnow().timetuple())

@register("gen/timestamp_str", str)
def timestamp_str():
  return datetime.now().ctime()

@register("build/dir", str)
def build_cwd():
  return os.getcwd()

@register("build/timestamp_str", rawstr)
def build_timestamp_str():
  return "__DATE__ \" \" __TIME__"

@register("build/hostname", str)
def build_hostname():
  return socket.gethostname()

@register("build/env", str)
def build_environment():
  out = []
  for key, val in os.environ.items():
    out.append(str(key) + "=" + str(val))
  return "\n".join(out)

@register("build/compiler_id", str)
def build_compiler():
  return os.environ["CMAKE_CXX_COMPILER_ID"]

@register("build/cmake_cache", str)
def build_cmake_cache():
  with open(os.environ['CMAKE_CACHEFILE'], 'r') as f:
    return "\n".join(f.readlines())

@register("build/built_by", str)
def build_built_by():
  return getpass.getuser()

@register("git/last_commit_hash", str)
def git_last_commit_hash():
  return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('utf8').strip()

@register("git/uncommitted_changes", bool)
def git_uncommitted_changes():
  return (subprocess.call(['git', 'diff', '--quiet']) != 0)

@register("git/diff_stat", str)
def git_diff_stat():
  lines = subprocess.check_output(['git', 'diff', '--stat']).decode('utf8').strip().split('\n')
  return '\n'.join(lines[:-1])  # chop off last line (summary)

@register("git/patch", str)
def git_patch():
  return subprocess.check_output(['git', 'diff-files', '-p']).decode('utf8')


#### -------------------------------------------------------------------------

# helper functions for output
def py_type_to_cpp_type(type):
  if type == str or type == rawstr:
    return 'const char*'
  elif type == uint64_t:
    return 'uint64_t'
  elif type == int:
    return 'int'
  elif type == bool:
    return 'bool'
  else:
    raise Exception("Unknown python-cpp type mapping for type " + type.__name__)

def sanitize_cpp_symbol(name):
  invalid_chars = '!@#$%^&*()-+=/\\[]{}`~.,;:\'"'
  for c in invalid_chars:
    name = name.replace(c, '')
  name = name.replace(' ', '_')

  if len(name) < 1:
    raise Exception("Invalid cpp symbol (contains no valid characters)")
  return name

def do_line_breaks(s, max_length):
  s = s.split('\n')

  # split into lines every 80 chars if possible
  split_after = [' ', ':', ';', ',']

  i = 0
  out = []
  while i < len(s):
    line = s[i]
    if len(line) < max_length:
      out.append(line)
      i += 1
    else:
      split_at = None
      for delim in split_after:
        idx = line.rfind(delim, 0, max_length)
        if idx > 0 and (split_at == None or idx > split_at):
          split_at = idx

      if split_at != None:
        out.append(line[0:split_at+1])
        s[i] = line[split_at+1:]
        continue
      else:
        # no good place to split, give up
        out.append(line)
        i += 1

  return out

def sanitize_cpp_value(type, val):
  if type == str:
    s = str(val)
    s = s.replace('\\', '\\\\')  # escape escapes
    s = s.replace('"', '\\"')  # escape quotes
    s = s.replace('\n', '\\n\n')  # escape newlines
    s = s.replace('\r', '\\r')  # escape newlines
    s = do_line_breaks(s, 80 - 6)  # s is now a list of lines

    out = '' if len(s[0]) < 50 else '\n    '
    out += '"' + '"\n    "'.join(s) + '"'
    return out

  elif type == rawstr:
    return str(val)
  elif type == int or type == uint64_t:
    return str(val)
  elif type == bool:
    return 'true' if val else 'false'
  else:
    raise Exception("Unknown cpp value sanitization for type " + type.__name__)

def cpp_missing_value(type):
  if type == str or type == rawstr:
    return 'NULL'
  elif type == int or type == uint64_t:
    return '0'
  elif type == bool:
    return 'false'
  else:
    raise Exception("Unknown cpp default value for type " + type.__name__)

### ---------------------------

parser = argparse.ArgumentParser(description="Generate build info for a C++ project.")
parser.add_argument('--disable', nargs='*', default=[], type=str, help="List of sections to not include.")
parser.add_argument('--skip', nargs='*', default=[], type=str, help="List of sections to still include, but use 'missing' values.")
parser.add_argument('--output', dest='output_file', type=argparse.FileType('w'), default='build_info.h', help="File to output to.")
parser.add_argument('--struct-name', type=str, default='BuildInfo', help="Name of the struct to generate. Should be a valid C++ symbol name.")
parser.add_argument('--helpers', type=bool, default=True, help="Include helper functions (like print).")

args = parser.parse_args()

# Only attempt git stuff if this folder has a git repo
def is_git_repo():
  return (subprocess.call(['git', 'rev-parse'], stdout=DEVNULL, stderr=subprocess.STDOUT) == 0)

if not is_git_repo():
  args.skip.append("git")

def disabled(category, entry):
  return category in args.disable or (category + '/' + entry) in args.disable

def skip(category, entry):
  return category in args.skip or (category + '/' + entry) in args.skip


### ----------------------------------

RED_COLOR = '\033[91m'
EMPHASIS_COLOR = '\033[1m'
END_COLOR = '\033[0m'

args.output_file.write('\n'.join([
  '#pragma once',
  '',
  '#include <string.h>',  # for strcmp
  '#include <ostream>',  # for std::ostream
  '',
  '// This is an automatically generated file.',
  '',
  'struct ' + args.struct_name + ' {',
  ''
  '  static constexpr const char* name = "' + args.struct_name + '";',
  ''
]))

missing_list = []
for category_name, category_val in categories.items():
  for entry_name, entry_val in category_val.items():
    if disabled(category_name, entry_name):
      continue

    type = entry_val[0]
    func = entry_val[1]

    cpp_type = py_type_to_cpp_type(type)
    cpp_name = sanitize_cpp_symbol(category_name + "_" + entry_name)
    cpp_value = None

    if not skip(category_name, entry_name):
      try:
        cpp_value = sanitize_cpp_value(type, func())
      except:
        sys.stderr.write("Unable to get value for " + RED_COLOR + category_name + "/" + entry_name + END_COLOR + ":\n")
        traceback.print_exc(file=sys.stderr)

    if cpp_value == None:
      cpp_value = cpp_missing_value(type)
      missing_list.append(category_name + "/" + entry_name)

    args.output_file.write('  static constexpr ' + cpp_type + ' ' + cpp_name + ' = ' + cpp_value + ';\n')

  # newline between sections
  args.output_file.write('\n')

# Helper functions
if args.helpers:
  # sections list
  args.output_file.write('  static constexpr const char* sections[] = {\n')
  for category_name, category_val in categories.items():
    for entry_name, entry_val in category_val.items():
      if disabled(category_name, entry_name):
        continue
      args.output_file.write('    "' + category_name + '/' + entry_name + '",\n')
  args.output_file.write('  };\n\n')

  # print by name
  args.output_file.write('  static void print(std::ostream& stream, const char* section = NULL) {\n')
  for category_name, category_val in categories.items():
    for entry_name, entry_val in category_val.items():
      if disabled(category_name, entry_name):
        continue

      cpp_name = sanitize_cpp_symbol(category_name + "_" + entry_name)
      name = category_name + '/' + entry_name

      args.output_file.write('    if (section == NULL || strcmp(section, "' + name + '") == 0)\n')
      if name in missing_list:
        args.output_file.write('      stream << (section == NULL ? "' + EMPHASIS_COLOR + name + END_COLOR + ': " : "") << "MISSING" << (section == NULL ? "\\n" : "");\n')
      else:
        args.output_file.write('      stream << (section == NULL ? "' + EMPHASIS_COLOR + name + END_COLOR + ': " : "") << ' + cpp_name + ' << (section == NULL ? "\\n" : "");\n')
  args.output_file.write('  }\n\n')

# close struct
args.output_file.write('};')

args.output_file.write('\n\n#ifdef ' + args.struct_name.upper() + '_IMPL\n')
args.output_file.write('constexpr const char* const ' + args.struct_name + '::sections[];\n')
args.output_file.write('#endif')
