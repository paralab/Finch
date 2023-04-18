TalyFEM 7.1.0 Release Notes
===========================

## New features

* Updated from PETSc 3.5 to PETSc 3.6.
* 64-bit PETSc indices work again.
* Added `ELEM::GetMeasure()` - the "measure" is the n-dimensional volume (i.e. length for 1D, area for 2D, and volume for 3D). This is currently only implemented for 1D and tetrahedral elements. (David)
* Added ElementTest tests.
  - These test `ELEM::Measure()`, `ELEM::CalculateNormal()`, and `ELEM::IsInnerPoint()`
* Added Basis tests.
  - Tests 1D/2D/3D linear/quadratic/cubic box, linear triangle, and linear tetrahedron. Thanks to Ramin and Aditya for the numbers.
* Enabled ParMETIS options to set a constant RNG seed. This should help reproducibility.
* Added `GRID::Scale` (David)
* Basis functions (and FEMElm) have been completely rewritten.
  - Proper cubic basis functions have been added for 1D/2D/3D box.
  - Triangle now uses an isoparametric mapping.
  - The new basis functions have been tested in two ways:
    + We've checked the error convergence for a steady-state heat problem for each basis function. This was done by plotting log(l2_error) vs log(n_elms) and looking at the slopes.
    + We've compared results with independent MATLAB codes developed by other people in the group (these make up the "BasisTest" test).
  - Normals are no longer calculated in the basis function (see below).
  - See "Breaking Changes" for changes to FEMElm.
* For box elements, higher-dimension volume and surface gauss points are automatically calculated from the 1D gauss points. You now only need to add gauss points/weights entry for 1D - 2D and 3D volumes/surfaces will be automatically defined for you.
* Triangle element connectivity order must be counter-clockwise. There is now a check for this. If an element isn't counter-clockwise, an exception will be thrown during loading.
* Tecplot reader now accepts the "LINE" element type for 1D. (Songzhe)
* Added support for checking test code coverage using lcov. (David)
  -  0) install lcov
  -  1) use gnu compiler
  -  2) run: `make coverage`
  -  3) results will be shown in taly_coverage/index.html
* Surface normals are now calculated with `ZEROPTV ELEM::CalculateNormal(const GRID* grid, int surface_id)`.
* Surface normals are now cached on the SurfaceIndicator object. Normals are calculated once, typically when the element is created. If you modify node positions and care about normals, be sure to call `ELEM::CalculateNormals(grid)` on any affected elements.
* Support for loading GMSH `*.msh` files.
  - Supports order 1/2/3 2D/3D box elements, as well as order 1 triangle and tetrahedron.
  - Supports domain decomposition.
  - Supprts loading node indicators from the Gmsh file, and uses the surface data from the Gmsh file when possible.
  - Does not support loading periodic data or element tags.
* Added `CEquation::Integrands4sideVectorOnly(const FEMElm& fe, int sideInd, ZEROARRAY<double>& be)` (equivalent of IntegrandsVectorOnly, but for surface integration)
* Added a KD-tree to `GRID` to speed up finding what elements are near a particular point. (Based on Mike Davies' code)
  - This should drastically speed up code using `GRID::FindElmAndLocPt`/`GridField::ValueAtPoint`.
  - If you want to directly query the tree yourself, use `GRID::kd_tree()` to get a reference to it. See `Grid/kdtree.h` for what kinds of queries you can do.
  - The tree is only constructed when it is first queried (i.e. when you first call `GRID::FindElmAndLocPt`). If you modify node points after causing the tree to be built, you should call `GRID::kd_tre().rebuild()` to rebuild the tree.
* `GRID::FindElmAndLocPt`/`GridField::ValueAtPoint` now use a two-phase system for finding the element containing a point, which should be much faster:
  - Broad-phase: query the KD-tree using the max radius of all elements to find elements near the point.
  - Narrow-phase: call `ELEM::IsInnerPoint` to see if the point is inside any element.
  - If no element is found, an exception is thrown.
* Added an implementation for `ELEM::IsInnerPoint` that should work for any **convex** element. If you have a non-convex element, you may get false negatives (in which case, an exception will be thrown).
* Config files are now read with [libconfig](http://www.hyperrealm.com/libconfig/) instead of Jaz.
* InputData boolean fields are now actually bools instead of ints.
* Added a lot of Doxygen documentation.
* Added reproducibility measures.
  - You can now use the "Repro" class (declared in Utils/reproducibility.h) to generate a 'repro.cfg' file for your project.

  - Basic usage:

    ```cpp
    int main(int argc, char** argv) {
      PetscInitialize(&argc, &argv, NULL, NULL);
      Repro r(argc, argv);
      r.write();  // must be AFTER PetscInitialize
    }
    ```

  - This will generate a repro.cfg file in the current working directory at program startup.
  The repro.cfg file contains things like TalyFEM build information (last git commit, git diff, compiler flags, linker flags, build directory, environment, ...) and run-time information (number of processes, config.txt, hostname, environment, command line arguments, ...).
  This information can be used to build a Pretty Good (TM) reconstruction of the TalyFEM version your code is using and the run environment, including minor uncommitted changes.

  - **IMPORTANT NOTE:** This will not capture any input files other than config.txt by default. To include other input files, you can use:

  ```cpp
  r.add_input_file("path.plt");
  ```
  
  This will include "path.plt" in the repro.cfg. The file is encoded in base64, so this can also be used for non-text (binary) files.

  - Advanced users can use the PythonScripts/build_info_gen.py script to generate their own BuildInfo struct and register it like so:

    ```cpp
    #define MYBUILDINFO_IMPL
    #include <MyBuildInfo.h>

    // ... before r.write() ...
    r.add_build_info<MyBuildInfo>();
    ```

    This will extend the source code reproducibility information to your own code, in addition to TalyFEM's.

## Breaking changes

* FEMElm has changed.
  - New constructor:
  
    ```cpp
    FEMElm fe(const GRID* grid, unsigned int flags = BASIS_ALL);
    ```

  - refill:

    ```cpp
    // OLD
    fe.refill(const GRID* pGrid, int elmID, const SurfaceIndicator* surface);

    // New (volume integration)
    fe.refill(const ELEM* elem, kBasisFunction bf, int relative_order_of_integration);
    fe.refill(int elemID, kBasisFunction bf, int rel_order);

    // New (surface integration)
    refill_surface(const ELEM* elem, const SurfaceIndicator* surface, kBasisFunction bf, int rel_order);
    refill_surface(int elemID, const SurfaceIndicator* surface, kBasisFunction bf, int rel_order);
    ```

    `bf` is a constant representing the type of basis function to use (`BASIS_LINEAR`, `BASIS_QUADRATIC`, ...) - you should probably pull this from InputData. `rel_order` affects how many gauss points to use. Valid values are -1/0/1, where -1 will use fewer gauss points, 0 will use the recommended number of gauss points, and 1 will use extra gauss points.

    The only thing that changes for surface integration is the gauss points.

  - The integration loop

    ```cpp
    // --------------------------------
    //   OLD
    // --------------------------------
    FEMElm fe;

    // for each element
    fe.refill(p_grid, element_id);
    fe.setRelativeOrder(order);)
    fe.initNumItg();
    while (fe.moreItgPts()) {
      fe.update4nextItgPt();
      // do calculations with fe
    }

    // --------------------------------
    //   NEW
    // --------------------------------
    FEMElm fe(p_grid);

    // for each element
    fe.refill(element_id, basis_function, rel_order);
    while (fe.next_itg_pt()) {
      // do calculations with fe
    }
    ```

  - Calculating basis functions at a user-defined integration point:

    ```cpp
    // OLD
    fe.setLocalEvalPt(loc_pt);
    fe.update4nextItgPt();
    // do calculations with fe

    // NEW
    fe.calc_at(ZEROPTV my_pt);
    // do calculations with fe
    ```

  - Important changed accessors:
    + `ZEROARRAY<double> fe.position` -> `ZEROPTV fe.position()` (integration point in global space)
    + `double detJxW` -> `double detJxW()` (weighted absolute value of jacobian determinant)
    + `double detJ()` -> `double jacc()` (unweighted jacobian determinant)
    + `const ELEM* fe.pElm` -> `const ELEM* fe.elem()`

  - Important new accessors:
    + `int fe.nbf()` - number of basis functions, use this in your loops instead of ELEM::n_nodes()!

  - Read Grid/femelm.h for more information. It is well-commented.

* CEquation::redimSolver has changed.
  - Now takes basis function and rel_order parameters:
  
    ```cpp
    // OLD
    equation.redimSolver(GRID* p_grid, int n_dof, const PeriodicData* periodic_data = NULL);

    // NEW
    equation.redimSolver(GRID* p_grid, int n_dof, kBasisFunction bf, int rel_order, const PeriodicData* = NULL);
    ```

  - Pull bf and rel_order from InputData:
    + `kBasisFunction bf = input_data.basisFunction;`
    + `int rel_order = input_data.basisRelativeOrder;`

  - This change was made to support hermite basis functions. Hermite basis functions have extra degrees of freedom at each node, which we need to know at redimSolver time.

* The order of gauss points has changed. This should only affect you if your code depends on the current integration point ID to do something.

* The order of some surface node points has changed. Surface nodes are now consistently counter-clockwise. The nodes in each surface have not changed, only their order has changed. Any code that uses `ELEM::GetSurfaceCheckArray()` or relies on a hard-coded order of surface node points may be affected by this.

* "orderOfBF" has been removed from InputData and config.txt.
  - It has been replaced with two new options:
  ```
  basisFunction = linear/quadratic/cubic/hermite (default: linear)
  basisRelativeOrder = -1/0/1 (default: 0)
  ```
  - `basisFunction` is which basis function to use. Not all basis functions have been implemented for all element types - an exception will be thrown if an option has not been implemented.
  - `basisRelativeOrder` is the relative order of integration. This only affects the number of gauss points.

* The `TezduyarUpwindFE` stabilizer terms are now zero-indexed.
  - SUPG(i) -> SUPG(i-1)
  - PSPG(i, j) -> PSPG(i-1, j-1)
  - DAPG(i) -> DAPG(i-1)

* Fixed off-by-one error in GRID::ProcIDOfElm.
  - GRID::ProcIDOfElm() ==> GRID::GetGridIdOfElementOwner

* ESSCondition::selfIntroduction() has been replaced by the stream operator (<<). Instead of:

  ```cpp
  ESSCondition boundary;
  boundary.selfIntroduction();
  ```
  use:
  
  ```cpp
  ESSCondition boundary;
  std::cout << boundary;
  ```

* ESSCondition::Specify() has been removed. Values are now set via the constructor.

    Instead of:
    
    ```cpp
    ESSCondition bc();
    bc.Apply(index, matrix_value, value);
    ```
    use:
    
    ```cpp
    ESSCondition bc(index, matrix_value, value);
    ```

* CEquation::initEssBC() no longer removes periodic bounds. If your code was relying on this, call CEquation::initPerBC() after calling initEssBC() to mimic the old behavior.

* If you modify node points after calling `GRID::FindElmAndLocPt` or `GridField::ValueAtPoint`, you must call `GRID::kd_tre().rebuild()` to rebuild the KD-tree - otherwise you will get incorrect results the next time you query the tree.

* The config file format has changed slightly since we switched to libconfig.
  - Booleans are now `True`/`False` instead of `1`/`0`.
  - Strings must be quoted.
  - If your code follows the typical input format used in the library, there should be four changes required:

    1) When reading values, do not pass the configuration mapping (which is no longer used):
    
    		ReadValue(conf, key, value) ==> ReadValue(key, value)

    2) If you override InputData::Initialize, change to the new signature:
    
    		Initialize(conf, cf) ==> Initialize()

    3) If you override ```ReadFromFile```, you should have code similar to:
    
    ```cpp
    // declare the map that stores key and values
    std::map < std::string, std::string > conf;
    jaz::ConfigFile cf;

    // read the input file using the jaz library
    if (cf.read(filename, conf) == false) {
      if (GetMPIRank() == 0)  // only rank 0 should print errors
        std::cerr << "Config.txt line: " << cf.error() << std::endl;
      return false;
    }

    // Initialize the basic fields
    InputData::Initialize(conf, cf);
    ```

    Replace all of this with:
    
    ```cpp
    // Read config file and initialize basic fields
    InputData::ReadFromFile(filename);
    ```

    If your code is similar to the above, except that you do NOT call
    
    `InputData::Initialize(conf, cf);`, instead replace it with:
    `InputData::ReadConfigFile(filename);`

    4) Update your config files:
    
    ```
    string input values MUST be surrounded by quotes:
    filename = out.plt  // NO!
    filename = "out.plt"  // YES
    ```

    bool input values MUST be changed to `True` or `False`:
    
    ```
    ifBoxGrid = 1  // NO!
    ifBoxGrid = True  // YES
    ```

  - The following items in InputData are now bool (previously int):
  
    ```cpp
    bool ifBoxGrid;
    bool ifTriElem;
    bool ifDD;
    bool ifPrintStat;
    bool ifPrintLog;
    bool ifPrintWarn;
    bool ifPrintInfo;
    bool ifPrintTime;
    bool ifLoadNodeIndicators;
    bool ifWriteNodeIndicators;
    ```

    This requires minor changes to user code:
    1) Tests of the values above need to change to reflect the new type. For example:
    
    ```cpp
    if (ifDD == 0) { ... }   ==>   if (!idDD) { ... }
    if (ifDD == 1) { ... }   ==>   if (idDD) { ... }
    ```

    2) The variables above were always assigned 1 or 0 in the library. If user code gave them different values for used the actual numeric values for some reason, that code will need to change.

    3) Configure files need to change to read in True or False (no quotes) instead of 1 or 0. For example:
    
    ```
    ifBoxGrid = 0  ==>  ifBoxGrid = False
    ifBoxGrid = 1  ==>  ifBoxGrid = True
    ```


## Bug fixes/misc

* Fixed many integer overflow bugs with large numbers of nodes (64-bit indices).
* Fixed get_nodes_in_element(kElem1d) ignoring basis function.
* Removed `xmls` from XDMF files, which fixes a ParaView warning when using HDF5.
* Many side effects of rewriting the basis functions...
  - Order "3" is now properly cubic, instead of quadratic as it was in the old library.
  - Tetrahedron surface gauss points have been corrected.
  - Surface normals should be consistently outward-facing now.
* GRID::operator= should copy all of its fields now. (Songzhe)
* All of the "find element containing point" methods in GRID/GridField
* Fixed "new nonzero" error for triangle elements (was due to a poor serial preallocator approximation). This may have a side effect of slightly increasing the preallocation for box elements.
* Serial and non-serial (no DD) now use the same code path in matrix creation/preallocation.


Building
========

## Dependencies

### PETSc 3.6.x

```bash
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.6.4.tar.gz
tar xvf petsc-3.6.4.tar.gz
rm petsc-3.6.4.tar.gz

cd petsc-3.6.4
export PETSC_DIR=`pwd`
./configure --with-debugging=1 --download-mpich --download-scalapack --download-parmetis --download-fblaslapack --download-blacs --download-metis --download-mumps
make

# Save PETSC_DIR
echo "PETSC_DIR=$PETSC_DIR" >> ~/.bashrc
```

### libconfig

TalyFEM 7.1.0 has one new dependency: [libconfig](http://www.hyperrealm.com/libconfig/). This is used for parsing config.txt files instead of the old jaz library.

**Installing libconfig from package manager (Ubuntu/Debian)**

If you are on a local Ubuntu install, you can just install libconfig as a package.

```bash
sudo apt-get install libconfig++-dev

echo "LIBCONFIG_DIR=/usr" >> ~/.bashrc
source ~/.bashrc
```

**Installing libconfig from source (Clusters)**

```bash
wget http://www.hyperrealm.com/libconfig/libconfig-1.5.tar.gz
tar xvf libconfig-1.5.tar.gz
rm libconfig-1.5.tar.gz

cd libconfig-1.5
export LIBCONFIG_DIR=`pwd`/install
./configure --prefix=$LIBCONFIG_DIR
make
make install

# Save LIBCONFIG_DIR and add the shared object file to LD_LIBRARY_PATH
echo "LIBCONFIG_DIR=$LIBCONFIG_DIR" >> ~/.bashrc
echo "LD_LIBRARY_PATH=\$LD_LIBRARYPATH:$LIBCONFIG_DIR/lib" >> ~/.bashrc

# Load update LD_LIBRARY_PATH to include libconfig++.so's directory
source ~/.bashrc
```


## TalyFEM

Set `PETSC_DIR=/path/to/petsc/root` and `LIBCONFIG_DIR=/path/to/libconfig/install` (typically in your `.bashrc`).

Next, set up your makefile.config:

```bash
cp makefile.config.base makefile.config
```

Uncomment the set of lines that matches your system (should be a set containing CC, CC_FLAGS, and TALYFEM_RUN_COMMAND).

After that, just run `make` to build the library.


TalyFEM 7.1.0 On Clusters
=========================

## CyEnce

Use the shared PETSc 3.6 and libconfig in `/work1/baskargroup/shared/`.

**Recommended environment (~/.bashrc):**

```bash
export LIBCONFIG_DIR=/work1/baskargroup/shared/libconfig-1.5/install
export LD_LIBRARY_PATH=$LIBCONFIG_DIR/lib:$LD_LIBRARY_PATH
export PETSC_DIR=/work1/baskargroup/shared/petsc-3.6.4
```

## Blue Waters

Use the `PrgEnv-gnu` programming environment module. Load the `cray-petsc/3.6.3.0` module for PETSc. Compile libconfig from source.

**Recommended environment (~/.bashrc):**

```bash
module swap PrgEnv-cray PrgEnv-gnu
module load cray-petsc/3.6.3.0

# Compile libconfig and set LIBCONFIG_DIR as shown in the "Dependencies" section.
```

## Stampede

Use the `petsc/3.6` module for PETSc. Use the `impi` module for MPI instead of mvapich2.

(Note: the default MPI module (mvapich2) has caused strange problems for me when running programs. Switching to `impi` fixed this for me.)

**Recommended environment (~/.bashrc):**

```bash
module load impi
module load petsc/3.6

# Compile libconfig and set LIBCONFIG_DIR as shown in the "Dependencies" section.
```

## Comet (SDSC)

The default compiler is Intel 2013. You need at least 2015 to build TalyFEM.

The provided PETSc 3.6 module was compiled with the 2013 edition of the Intel compiler. TalyFEM must be compiled with 2015 or later (due to the use of C++11 code). If you try to use the system PETSc, you will run into linker errors: `undefined reference to '__intel_avx_rep_memset'`.

This happens because the PETSc library flags are causing the compiler to try and link with the 2013 version of the Intel runtime (libirc), which does not have `__intel_avx_rep_memset`.

The simplest work-around for this is to compile PETSc 3.6 from source yourself. You can use the system MPICH with PETSc by loading mvapich2 (`module load mvapich2_ib`) and configuring PETSc with `--with-mpi=/opt/mvapich2/intel/ib` instead of `--download-mpich`.

**Recommended environment (~/.bashrc):**

```bash
module swap intel intel/2016.3.210
module load mvapich2_ib

# Compile PETSc and set PETSC_DIR as shown in the "Dependencies" section.
# compile libconfig and set LIBCONFIG_DIR as shown in the "Dependencies" section.
```
