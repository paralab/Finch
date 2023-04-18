#!/usr/bin/env python

import subprocess
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import libconf
import sys

if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
  print("Invalid path to SSHT executable (argument 1).")
  sys.exit(1)

PROGRAM_PATH = sys.argv[1]  # path to SSHT executable

def write_mesh(nsd, order_num, n_elems, mesh_path, output_path):
  with open('/dev/null', 'w') as devnull:
    cmd = ["gmsh", mesh_path, "-o", output_path, "-" + str(nsd), "-order", str(order_num), "-setnumber", "cl", str(1.0 / n_elems)]
    print(cmd)
    subprocess.check_call(cmd, stdout=devnull)

def get_error(cfg):
  # write config.txt using the libconf library
  with open('config.txt', 'w') as f:
    libconf.dump(cfg, f)

  # run the program with the direct solver
  #out = subprocess.check_output(['mpirun', '-n', '8', PROGRAM_PATH])
  #out = subprocess.check_output(['mpirun', '-n', '7', '-bind-to', 'hwthread', PROGRAM_PATH, '-pc_type', 'lu', '-pc_factor_mat_solver_package', 'mumps', '-ksp_type preonly'])
  out = subprocess.check_output(['mpirun', '-n', '7', '-bind-to', 'hwthread', PROGRAM_PATH, '-ksp_rtol', '1e-10', '-ksp_atol', '1e-10'])

  # parse the error using a regular expression
  match = re.search("L2 error = ([\d.e-]+)", out)
  err = float(match.group(1))

  # make sure the error isn't 0 - it should never be this low.
  # if it is, something went wrong somewhere
  if err < 1e-16:
    raise RuntimeError("Error too close to zero: " + str(err))

  return err

def plot_conv(cfgs, errors):
  # calculate our x and y values for the plot/linear regression
  # (make sure to use floating point divide, not integer divide)
  x = [math.log10(cfg['Nelemx']) for cfg in cfgs]
  y = [math.log10(e) for e in errors]

  # use numpy to fit a line so we can see the slope (m = slope, b = y-intercept)
  m, b = np.polyfit(x, y, 1)
  plt.plot(x, y, marker='.')
  plt.xlabel("log10(Lx/Nelemx)")
  plt.ylabel("log10(l2_error)")
  plt.title("Slope: {0:.8f}".format(m))
  plt.show()

def plot_all(results):
  labels = []
  for name, result in results.iteritems():
    x = [math.log10(n) for n in result[0]]
    y = [math.log10(e) for e in result[1]]
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, y, marker='.')
    labels.append((name + ' ({0:.4f})').format(m))

  plt.title("TalyFEM Convergence")
  plt.xlabel("log10(Nelemx)")
  plt.ylabel("log10(l2_error)")
  plt.legend(labels)
  plt.show()

if __name__ == "__main__":
  # Main program
  neumann_lr = {"boundaries": {"left": "neumann", "right": "neumann"}}
  neumann_tb = {"boundaries": {"top": "neumann", "bottom": "neumann"}}
  tests = {
    "linear triangle": (2, "linear", [10, 20, 40, 80], {}),
    "linear triangle (neumann l/r)": (2, "linear", [10, 20, 40, 80], neumann_lr),
    "quadratic triangle": (2, "quadratic", [10, 20, 40, 80], {}),
    "quadratic triangle (neumann l/r)": (2, "quadratic", [10, 20, 40, 80], neumann_lr),
    "quadratic triangle (neumann t/b)": (2, "quadratic", [10, 20, 40], neumann_tb),
    "cubic triangle": (2, "cubic", [10, 20, 40, 80], {}),
    "cubic triangle (neumann l/r)": (2, "cubic", [10, 20, 40, 80], neumann_lr),
    "cubic triangle (neumann t/b)": (2, "cubic", [10, 20, 40], neumann_tb),

    "linear tetrahedron": (3, "linear", [10, 20, 40], {}),
    "linear tetrahedron (neumann l/r)": (3, "linear", [10, 20, 40], neumann_lr),

    "quadratic tetrahedron": (3, "quadratic", [10, 20, 40], {}),
    "quadratic tetrahedron (neumann l/r)": (3, "quadratic", [10, 20, 40], neumann_lr),

    #"cubic tetrahedron": (3, "cubic", [6, 12, 24], {}),
    #"cubic tetrahedron (neumann l/r)": (3, "cubic", [6, 12, 24], neumann_lr),
  }

  results = {}
  for testname, test in tests.iteritems():
    nsd = test[0]
    basis = test[1]
    n_elems = test[2]

    if basis == "linear":
      order_num = 1
    elif basis == "quadratic":
      order_num = 2
    elif basis == "cubic":
      order_num = 3
    else:
      raise Exception("Unknown basis")

    cfgs = []
    errors = []
    results[testname] = ([], [])
    print testname
    for n in n_elems:
      base_msh = "tri_base.geo"
      if nsd == 3:
        base_msh = "tet_base.geo"

      if os.path.exists("mesh.msh"):
        os.remove("mesh.msh")  # remove old file in case write_mesh fails, so we'll notice

      write_mesh(nsd, order_num, n, base_msh, "mesh.msh")

      # describe the contents of config.txt as a Python dictionary
      cfg = {
        "ifDD": True,
        "nsd": nsd,
        "basisFunction": basis,
        "ifBoxGrid": False,
        "inputFilenameGrid": "mesh.msh",
        "ifLoadNodeIndicators": True,
        "Lx": 1.0,
        "Nelemx": n,
      }

      cfg.update(test[3])

      err = get_error(cfg)
      cfgs.append(cfg)
      errors.append(err)
      results[testname][0].append(n)
      results[testname][1].append(err)

    # show the results
    print(testname, "errors: " + ", ".join(["{0:.6E}".format(e) for e in errors]))
    # plot_conv(cfgs, errors)

  plot_all(results)
