#!/usr/bin/env python

import os
import sys
import itertools  # for product

sys.path.append(os.environ['TALYFEM_DIR'] + '/python_scripts')

from test_common import run_tests
from test_common import DEFAULT_PROGRAM_PATH
from test_common import DEFAULT_PROGRAM_ARGS
from test_common import have_mumps
from test_common import all_tests_passed

# --- base tests ---

generated_mesh_tests = [
{
  # test_0
  "__DESCR": "1D, generated mesh",
  "nsd": 1,
  "ifBoxGrid": True,
  "Lx": 1,
  "Nelemx": 48,
}, 
{
  # test_1
  "__DESCR": "2D, generated mesh",
  "nsd": 2,
  "ifBoxGrid": True,
  "Lx": 1,
  "Ly": 1,
  "Nelemx": 24,
  "Nelemy": 24,
}, 
{
  # test_2
  "__DESCR": "3D, generated mesh",
  "nsd": 3,
  "ifBoxGrid": True,
  "Lx": 1,
  "Ly": 1,
  "Lz": 1,
  "Nelemx": 24,
  "Nelemy": 24,
  "Nelemz": 24,
}, 
{
  # test_3 (triangles)
  "__DESCR": "2D tri, generated mesh",
  "nsd": 2,
  "ifBoxGrid": True,
  "ifTriElem": True,
  "Lx": 1,
  "Ly": 1,
  "Nelemx": 24,
  "Nelemy": 24,
}, 
]

gmsh_meshes = [
  ("1D_linear_x.msh", "linear", "line_x"),
  ("1D_linear_y.msh", "linear", "line_y"),
  ("1D_linear_z.msh", "linear", "line_z"),

  ("1D_quad_x.msh", "quadratic", "line_x"),
  ("1D_quad_y.msh", "quadratic", "line_y"),
  ("1D_quad_z.msh", "quadratic", "line_z"),

  ("1D_cubic_x.msh", "cubic", "line_x"),
  ("1D_cubic_y.msh", "cubic", "line_y"),
  ("1D_cubic_z.msh", "cubic", "line_z"),

  ("2D_linear_xy.msh", "linear", "plane_xy"),
  ("2D_linear_xz.msh", "linear", "plane_xz"),
  ("2D_linear_yz.msh", "linear", "plane_yz"),

  ("2D_quad_xy.msh", "quadratic", "plane_xy"),
  ("2D_quad_xz.msh", "quadratic", "plane_xz"),
  ("2D_quad_yz.msh", "quadratic", "plane_yz"),

  ("2D_cubic_xy.msh", "cubic", "plane_xy"),
  ("2D_cubic_xz.msh", "cubic", "plane_xz"),
  ("2D_cubic_yz.msh", "cubic", "plane_yz"),

  ("2D_linear_xy_tri.msh", "linear", "plane_xy"),
  ("2D_linear_xz_tri.msh", "linear", "plane_xz"),
  ("2D_linear_yz_tri.msh", "linear", "plane_yz"),

  ("2D_sphere_octant.msh", "linear", "sphere"),
  ("cube_tet.msh", "linear", "cube"),
]

def build_gmsh_test(t):
  return {
    "__DESCR": "Gmsh mesh " + t[0],
    "__COPY_FILES": ["gmsh/" + t[0]],
    "__NUM_PROCESSORS": 1,
    "ifBoxGrid": False,
    "inputFilenameGrid": t[0],
    "basisFunction": t[1],
    "nsd": 3,
    "analyticSolution": t[2],
    "ifLoadNodeIndicators": True,
  }

gmsh_tests = map(build_gmsh_test, gmsh_meshes)

# --- modifiers ---

def test_basis_functions(base):
  if "ifTriElem" in base and base["ifTriElem"] == True:
    return [
      dict(base.items() + {"basisFunction": "linear"}.items()),
    ]
  else:
    return [
      dict(base.items() + {"basisFunction": "linear"}.items()),
      dict(base.items() + {"basisFunction": "quadratic"}.items()),
      dict(base.items() + {"basisFunction": "cubic"}.items()),
    ]

def get_test_nsd(base):
  if "analyticSolution" in base:
    if base["analyticSolution"].startswith("line"):
      return 1
    elif base["analyticSolution"].startswith("plane"):
      return 2
    elif base["analyticSolution"].startswith("cube") or base["analyticSolution"].startswith("sphere"):
      return 3
    else:
      raise Exception("Unknown nsd")
  else:
    return base["nsd"]

# test *all possible* boundary conditions (combinatorial explosion!)
def test_boundaries_exhaustive(base):
  nsd = get_test_nsd(base)

  # exhaustive boundary tests
  conditions = ["dirichlet", "neumann"]
  boundaries = ["left", "right"]
  if nsd >= 2:
    boundaries += ["bottom", "top"]
  if nsd >= 3:
    boundaries += ["back", "front"]

  # octant
  if base.get("analyticSolution") == "sphere":
    boundaries = ["left", "bottom", "back"]

  # cartesian product of cartesian product...
  # *[list comprehension] unpacks the generated list into separate arguments to the product function
  additional = itertools.product(*[itertools.product([b], conditions) for b in boundaries])
  out = []
  for i in additional:
    has_dirichlet = False
    for boundary, condition in i:
      if condition == "dirichlet":
        has_dirichlet = True
        break

    if has_dirichlet:
      out.append(dict(base.items() + {"boundaries": dict(i)}.items()))
  return out

# test a much smaller subset of possible boundary condition combinations
def test_boundaries(base):
  nsd = get_test_nsd(base)

  additional = []
  if nsd >= 1:
    additional += [
      [], # all dirichlet
      [("left", "neumann"), ("right", "dirichlet")],
      [("left", "dirichlet"), ("right", "neumann")],
    ]

  if nsd >= 2:
    additional += [
      [("left", "dirichlet"), ("right", "neumann"), ("top", "neumann"), ("bottom", "neumann")],
      [("left", "neumann"), ("right", "neumann"), ("top", "dirichlet"), ("bottom", "neumann")],
    ]

  if nsd >= 3:
    additional += [
      [("left", "dirichlet"), ("right", "neumann"), ("top", "neumann"), ("bottom", "neumann"), ("front", "neumann"), ("back", "neumann")],
      [("left", "neumann"), ("right", "neumann"), ("top", "neumann"), ("bottom", "neumann"), ("front", "dirichlet"), ("back", "neumann")],
    ]

  out = []
  for i in additional:
    out.append(dict(base.items() + {"boundaries": dict(i)}.items()))
  return out

def test_parallel(base):
  nsd = get_test_nsd(base)
  if nsd >= 2 and base.get("ifTriElem", False) != True:
    return [
      dict(base.items() + {"ifDD": True, "__NUM_PROCESSORS": 4}.items()),   # parallel with DD
      dict(base.items() + {"ifDD": False, "__NUM_PROCESSORS": 4}.items()),  # parallel w/o DD
      dict(base.items() + {"ifDD": False, "__NUM_PROCESSORS": 1}.items()),  # serial
    ]
  else:
    # no DD in 1D or triangle case
    return [
      dict(base.items() + {"ifDD": False, "__NUM_PROCESSORS": 1}.items()),  # serial
    ]

def test_n_gps(base):
  return [
    dict(base.items() + {"basisRelativeOrder": 0}.items()),
    dict(base.items() + {"basisRelativeOrder": 1}.items()),
  ]

def map_and_flatten(func, lst):
  return [item for sublist in map(func, lst) for item in sublist]

# -- apply modifiers ---

generated_mesh_tests = map_and_flatten(test_boundaries_exhaustive, generated_mesh_tests)
generated_mesh_tests = map_and_flatten(test_basis_functions, generated_mesh_tests)
generated_mesh_tests = map_and_flatten(test_parallel, generated_mesh_tests)
#generated_mesh_tests = map_and_flatten(test_n_gps, generated_mesh_tests)

gmsh_tests = map_and_flatten(test_boundaries_exhaustive, gmsh_tests)
gmsh_tests = map_and_flatten(test_parallel, gmsh_tests)
#gmsh_tests = map_and_flatten(test_n_gps, gmsh_tests)

print len(generated_mesh_tests), "library-generated mesh tests"
print len(gmsh_tests), "gmsh-based tests"

PROGRAM_PATH = DEFAULT_PROGRAM_PATH
EXECUTABLE_PATH = os.environ['EXECUTABLE_PATH']
PROGRAM_ARGS = DEFAULT_PROGRAM_ARGS + ["-n", "__NUM_PROCESSORS", EXECUTABLE_PATH]

# skip mumps tests
#if not have_mumps():
#  for test in tests:
#    if "__PROGRAM_ARGS" in test and "mumps" in test["__PROGRAM_ARGS"]:
#      test["__SKIP"] = "MUMPS not installed"

tests = generated_mesh_tests + gmsh_tests
results = run_tests(tests, PROGRAM_PATH, PROGRAM_ARGS)
sys.exit(0 if all_tests_passed(results) else 1)
