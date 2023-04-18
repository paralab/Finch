# New Features

* Changed the build system to CMake - see [BUILDING.md](https://bitbucket.org/baskargroup/taly_fem/src/master/docs/BUILDING.md?at=master&fileviewer=file-view-default) for instructions on building TalyFEM, see [CMAKE.md](https://bitbucket.org/baskargroup/taly_fem/src/master/docs/CMAKE.md?at=master&fileviewer=file-view-default) for CMake tips and info on how to update your project (1cb116a48f14) (1424f048e6e38)
* Support for using [PTSCOTCH](https://www.labri.fr/perso/pelegrin/scotch/) instead of ParMETIS for domain decomposition (3f9acbb4335)
* Support for 1D and 2D elements defined in 3D space - use the `BASIS_DIMENSION_REDUCTION` basis flag with FEMElm/CEquation::add_basis_flag() to enable (af029ad9b211)
* Better support for 1D surface integration - 1D basis functions now use a "Point" basis function for surface integration that has only N and Jacobian values (b057dca2128d)
* Added `make docs` target that generates Doxygen documentation (694cd2a39b9)
* Automatically run `clang-tidy` during compilation if it is installed (cdf602008)
* Added quadratic and cubic (10-node) triangle basis functions (Pengfei) (note: surface integration may be incorrect for these elements)
* Added quadratic and cubic tetrahedron basis functions (Pengfei) (note: surface integration may be incorrect for these elements)
* `FEMElm::refill()` now infers the basis function from the element instead of requiring the user to specify it (but you can still do it the old way, e.g. to force linear) (e641adb6d1)
* Non-DD Gmsh grid loader now supports loading mixed element types in a mesh (e.g. a mesh with both triangles and quadrilaterals) (1e68c8f280b9d)
* Non-DD Tecplot output now supports writing meshes with mixed element types (by writing degenerate elements) (1e68c8f280b9d)
* Added support for relative order -2 and -3 for quadratic and cubic elements (3da5a914aec7)
* `GRID::GetLocalPtv` (the function that maps a physical point to an isoparametric point) now uses an iterative method. The returned point is guaranteed to have L2 error of `1e-7` or less. This fixes `GetLocalPptv` being incorrect with non-parallelogram elements. Based on Yu Xie's streamtracing code. (a7634bd699f5)
* Added experimental generic d2N calculation which should work for all basis functions - use the `BASIS_SECOND_DERIVATIVE` basis flag with FEMElm/CEquation::add_basis_flag() to enable (6961032bf8)
* Support for custom separators in the `Print*` functions (32043509257)
* Added `InputData::ReadArray`
* Added `FunctionIntegrator` class

# New Tests

* Add tests for `GRID::GetLocalPtv` for all element types (5e972a3ea5)
* Added steady-state heat tutorial which tests surface integration (Neumann BC) and out-of-plane elements (1D/2D elements in 3D space) for all element types (line/quad/tri/hex/tet) and orders (linear/quadratic/cubic) (38ee7019ac2)
* Added a test for the stabilizer (SUPG/PSPG) for line/quad/hex linear/quadratic/cubic and tri/tet linear elements/bases (61bed342e2e7b)
* Added TransientHeat tutorial cases for loading PLT files with Windows line endings (6aea9c4)
* Added all the tutorial/tests to CTest (used by `make test`) (77a9651a2f) (f66bfc0fcf36298)

# Breaking Changes

* Changed the build system to CMake - see [BUILDING.md](https://bitbucket.org/baskargroup/taly_fem/src/master/docs/BUILDING.md?at=master&fileviewer=file-view-default) for instructions on building TalyFEM, see [CMAKE.md](https://bitbucket.org/baskargroup/taly_fem/src/master/docs/CMAKE.md?at=master&fileviewer=file-view-default) for CMake tips and info on how to update your project (1cb116a48f14) (1424f048e6e38)
* Surfaces loaded from Gmsh files are now defined by Gmsh physical ID - to my knowledge, everyone already using Gmsh has already adjusted to this change (4dec85bc7121dd9b5)
* Require Python 2.7.x (for libconf.py) - you may need to load a Python module on some clusters now (ba386adee6fc49510)
* Removed `OVERRIDE` and `STATIC_ASSERT` macros now that we require C++11 - if you used these macros in your own code, you'll need to replace them with the real C++ override/static_assert keywords (d7c5594a2ea)

# Miscellaneous

* Surface integration now matches how it was done in 7.0 and previous, but the Jacobian is calculated using a rotated lower-dimension element (04c93bcdc7c11b9)
* Changed arguments for the `FEMElm::d2N` accessor from `(bf, i)` to `(bf, dir1, dir2)`. While `d2N` is logically a matrix, it has several repeated values, so it is stored in a "reduced" format that does not contain the duplicate values. This accessor abstracts that away, so `fe.d2N(i, ...)` acts like a matrix. (90658872a6c)
* Added bounds checking to FEMElm N/dN accessors (6eac38be837)
* 1D node ordering and basis functions now follow the same ordering conventions as 2D and 3D (first node is leftmost, second node is rightmost, 3+ fill in the center of the line from left to right) (f30d1fac99f)
* Replace default `CEquation::IntegrandsPreallocator` with Ae.fill(true), which is currently the only sane option - you shouldn't override IntegrandsPreallocator (9e4941a8e23)
* Added a few CLion project settings in the `.idea` directory for shared code style settings (51bb13764c0)
* Testing framework `test_common.py` now uses [libconf](https://pypi.python.org/pypi/libconf) to write config files in order to support more advanced config files (e.g. the new SSHT tutorial's "boundaries" structure). No need to install it, libconf is included in `python_scripts/external`. (e0c7af2c93a)
* Python scripts shebang line now uses `/usr/bin/env python` instead of `/usr/bin/python` to support loading python as a module (5b75cd83b7b46) (6c46e919d)
* The special `kTriangle2dIn3d` element type has been removed, since its functionality has been absorbed into the basis functions (a63cfeaed2229)
* Added generic d2Nde calculation to `NonlinearBasis` (097af1ec79)
* Properly escape carriage returns (\r) in build info (d01fe5c9df0)
* CEquation::redimSolver now takes a `bool hermite` parameter instead of a basis function, since the basis function can vary throughout the mesh now (hermite is not supported with mixed meshes). An overload that matches the old signature still exists. (e641adb6d1)
* Remove unnecessary `std::sort` calls in `KDTree::get_elms_near_pt`, which improves performance significantly (d26902b93)
* GridField no longer requires `NodeData::UpdateDataStructures()` to exist (3bcf90e3)
* Cache rotation matrix after first basis calculation
* Throw an exception if an invalid HDF5 dataset name is used as a NodeData name when writing HDF5 files (e89d569197b)
* Handle Windows line endings (\r\n) and tab characters when reading Tecplot headers (6aea9c447d)
* You can now change the basis flags that are passed to the FEMElm constructor inside CEquation::Assemble() with `CEquation::set_basis_flags()`, `add_basis_flag()`, and `remove_basis_flag()` (get with `CEquation::basis_flags()`)