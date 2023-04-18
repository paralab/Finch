# TalyFEM 7.3.0 Release Notes

## Known Issues

* Surface integration for quadratic/cubic tetrahedrons does not work correctly.
* Additionally, volume integration for cubic tetrahedrons may not work correctly. `fe.jacc()` returns a negative Jacobian determinant, however `fe.detJxW()` takes the absolute value, which appears to work correctly (according to SSHT convergence analysis). The cause is likely that the shape functions are associated with the wrong nodes.

## Changes
* **Minor breaking change:** PeriodicBounds object is no longer stored in CEquation. If you have any references to `ceqn->periodic_bounds_` in your code (you probably don't), change it to `ceqn->periodic_data_->periodic_bounds()`.
* `ceqn->PresetPeriodicData(...)` is now automatically called by CEquation (via redimSolver) and no longer needs to be called by user code. You may safely remove calls to `ceqn->PresetPeriodicData(...)` from your code.

## New features
* Added Mesquite support (by way of a Mesquite adapter that works with TalyFEM GRIDs). To use, pass `-DENABLE_MESQUITE=YES` to your CMake options. To specify the Mesquite search directory, you can pass `-DMESQUITE_DIR=/path/here` or set the `MESQUITE_DIR` environment variable (e.g. in your bashrc). See [docs/BUILDING.md](docs/BUILDING.md) for more info.
* Added support for higher-order triangles and tetrahedrons in domain decomposition.
* Some documentation (WIP): https://www.sharelatex.com/project/59d50d2e9a7bf21ed17eb923
* The HDF5 interface has been changed to support writing non-contiguous element data.

## Fixes 
* Fixed parallel-but-not-DD SSHT error calculation being multiplied by the number of processes.
* Fixed surface integration for higher-order triangles (quadratic/cubic), see below.
* Various fixes for higher-order triangles and tetrahedrons:
  - Fixed higher-order 2D triangle gauss point weights.
  - Added proper quadratic tetrahedron volume gauss points (8 instead of 10).
  - Fixed tetrahedron surface 3 (XZ-plane) surface gauss points (counter-clockwise X rotation should have been clockwise).
  - Fixed tetrahedron surface 4 (XY-plane) gauss point order (by mirroring along the XY axis, i.e. swapping Xs with Ys).
  - Fixed quadratic tetrahedron shape function order (previously shape function order did not match node order).
  - Note that surface integration for higher-order tetrahedrons still does not work correctly (see known issues).
