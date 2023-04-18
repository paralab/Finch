#!/usr/bin/python

import os
import sys
import subprocess
import re
import math
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd() + "/../../../PythonScripts")
from test_common import write_config

base_cfg = {
  "__NUM_PROCESSORS": 8,
  "ifDD": False,
  "typeOfIC": 1,
  "Lx": 1,
  "Ly": 1,
  "Lz": 1,
  "dt": 0.001,
  "nOfTS": 10,
}

cfgs = {
  "linear 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "linear" },
  "quadratic 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "quadratic" },
  "cubic 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "cubic" },

  "linear 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "linear" },
  "quadratic 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "quadratic" },
  "cubic 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "cubic" },

  "linear 2D triangle": { "nsd": 2, "ifBoxGrid": True, "ifTriElem": True, "basisFunction": "linear" },
}

n_elems = [12, 24, 48, 72]

all_tests = {}
for name, mod in cfgs.iteritems():
  base = base_cfg.copy()
  base.update(mod)

  test = []
  for n in n_elems:
    cfg = base.copy()
    cfg["Nelemx"] = n
    cfg["Nelemy"] = n
    cfg["Nelemz"] = n
    test.append(cfg)

  all_tests[name] = test


def get_error(cfg):
  cwd = os.getcwd()
  #os.chdir('/tmp/')
  write_config('config.txt', cfg)
  executable = cwd + "/../CODE/BT"
  out = subprocess.check_output(["mpirun", "-n", str(cfg["__NUM_PROCESSORS"]), executable, "-ksp_type", "preonly", "-pc_type", "lu", "-pc_factor_mat_solver_package", "mumps"])
  match = re.search("error = ([\d.e-]+)", out)
  #os.chdir(cwd)
  return float(match.group(1))

def plot(cfgs, errors):
  x = [math.log(cfg['Nelemx']) for cfg in cfgs]
  y = [math.log(e) for e in errors]
  m, b = np.polyfit(x, y, 1)
  plt.plot(x, y, marker='.')
  plt.xlabel("log(n_elems)")
  plt.ylabel("log(error)")
  plt.title("Slope: " + str(m))
  plt.show()

def plot_all(results):
  labels = []
  for name, result in results.iteritems():
    x = [math.log(n) for n in result[0]]
    y = [math.log(e) for e in result[1]]
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, y, marker='.')
    labels.append((name + ' ({0:.2f})').format(m))

  plt.title("TalyFEM Convergence")
  plt.xlabel("log(n_elems)")
  plt.ylabel("log(error)")
  plt.legend(labels)
  plt.show()

# Main program
results = {}
for name, cfgs in all_tests.iteritems():
  print name, "..."
  results[name] = ([cfg['Nelemx'] for cfg in cfgs], [get_error(cfg) for cfg in cfgs])
  print results[name], "\n"
  #plot(test, errors)

plot_all(results)

