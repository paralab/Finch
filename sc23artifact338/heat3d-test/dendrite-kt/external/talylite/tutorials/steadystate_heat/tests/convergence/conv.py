#!/usr/bin/python

import os
import sys
import subprocess
import re
import math
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd() + "/../../../../python_scripts")
from test_common import write_config

base_cfg = {
  "ifDD": False,
  "ifBoxGrid": True,
  "Lx": 1.0,
  "Ly": 1.0,
  "Lz": 1.0,
}

cfgs = {
  "linear 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "linear" },
  "quadratic 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "quadratic" },
  "cubic 2D box": { "nsd": 2, "ifBoxGrid": True, "basisFunction": "cubic" },
  "linear 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "linear" },
  "quadratic 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "quadratic" },
  "cubic 3D hexahedron": { "nsd": 3, "ifBoxGrid": True, "basisFunction": "cubic" },

  #"linear 2D triangle": { "nsd": 2, "ifBoxGrid": True, "ifTriElem": True, "basisFunction": "linear" },
}

n_elems = [6, 12, 24, 48, 72]

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
  executable = "/home/lofquist/talyfemlib/cmake-build-debug/tutorials/steadystate_heat/ssht"
  out = subprocess.check_output(["mpirun", "-n", '4', executable, "-ksp_rtol", '1e-10'])
  match = re.search("L2 error = ([\d.e-]+)", out)
  err = float(match.group(1))
  if err < 1e-16:
    raise Exception("Error too close to zero: " + str(err))
  return err

def plot(cfgs, errors):
  x = [math.log10(cfg['Lx'] / cfg['Nelemx']) for cfg in cfgs]
  y = [math.log10(e) for e in errors]
  m, b = np.polyfit(x, y, 1)
  plt.plot(x, y, marker='.')
  plt.xlabel("log(Lx/Nelemx)")
  plt.ylabel("log(l2_err)")
  plt.title("Slope: " + str(m))
  plt.show()

def plot_all(results):
  labels = []
  for name, result in results.iteritems():
    x = [math.log10(1.0 / n) for n in result[0]]
    y = [math.log10(e) for e in result[1]]
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, y, marker='.')
    labels.append((name + ' ({0:.4f})').format(m))

  plt.title("TalyFEM Convergence")
  plt.xlabel("log10(Lx/Nelemx)")
  plt.ylabel("log10(l2_error)")
  plt.legend(labels)
  plt.show()

# Main program
results = {}
for name, cfgs in all_tests.iteritems():
  results[name] = ([cfg['Nelemx'] for cfg in cfgs], [get_error(cfg) for cfg in cfgs])
  #plot(test, errors)

plot_all(results)

