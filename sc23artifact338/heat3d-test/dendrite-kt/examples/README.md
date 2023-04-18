Coding Conventions
===================
* The equation class and corresponding solver must start with same name. Example: `htEq` for equation class and `htSolver` is
correct whereas `talyEq` and `htSolver` combination is not. Correspondingly, the same prefix must be set as a solver option for PETSc.
* Make variable name as explicit as possible.
* Do not read dimensions from the file `idata->nsd` is not recommended, rather use `DIM` that is defined during compile time.
This will allow compiler to optimize the loop.

Regression tests
============
Mpi runner is by default as 
```text
runner=$(which mpirun)
```
Running tests for all 
```bash
sh regression_test.sh
```
Running tests for single example 
```bash
sh regression_test.sh {FOLDERNAME} # sh regression_test.sh TSHT 
```