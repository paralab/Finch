Every HPC environment has its own quirks. Here are some that our group has found.

## CyEnce (Iowa State University cluster)

### Compiling PETSc With Intel MPI Requires Specifying Compilers

Condo and CyEnce use non-standard names for the Intel MPI compiler wrapper (`mpiicc`, `mpiicpc`, and `mpiifort`). If you want to use the system MPI, you'll need to add `--with-cxx=mpiicpc`, `--with-cc=mpiicc`, and `--with-fc=mpif90` to your PETSc configure options.

#### 64-bit Indices on CyEnce/Condo/Stampede (Likely Out Of Date)

I kept getting weird Fortran errors about `mpi.mod` not being found when compiling. This seemed to fix it.

```bash
export I_MPI_FC=ifort
./configure --useThreads=1 --with-debugging=1 --with-cc=mpicc --with-cxx=mpicxx --with-clanguage=cxx --with-fc=mpifc --download-parmetis=1 --with-scalar-type=real --download-metis=1 --download-fblaslapack --with-64-bit-indices
```

`mpifc` is a script that runs `gfortran` by default, so we need to set `I_MPI_FC=ifort` to use the Intel Fortran compiler so that Fortran MPI works.


## Condo

### Compiling PETSc With Intel MPI Requires Specifying Compilers

Condo and CyEnce use non-standard names for the Intel MPI compiler wrapper (`mpiicc`, `mpiicpc`, and `mpiifort`). If you want to use the system MPI, you'll need to add `--with-cxx=mpiicpc`, `--with-cc=mpiicc`, and `--with-fc=mpiifort` to your PETSc configure options.

### Libconfig ./configure mkdir fails

On Condo, `$TMPDIR` is set to `/scratch/$USER`. This directory is not automatically created. Use `mkdir $TMPDIR` to create it, then re-run `./configure ...`.

### Libconfig Make has Linker Errors

When you run `make` for libconfig, it may give linker errors when using the Intel compiler on Condo. The cause seems to be that the libconfig example programs are failing to link with the Intel runtime - but the libconfig library itself still compiles fine. It is safe to ignore this error (and continue to run `make install`).

### Libconfig has strange YACC-related errors

Try configuring with `./configure YACC=no ...`.


## Blue Waters

#### Using Cray PETSc

Cray supplies a tuned version of PETSc for its systems. You may get better performance using the pre-compiled cray-petsc module instead of building it yourself.

```bash
# the module PrgEnv-cray is loaded by default
module load cray-petsc/3.6
```

#### Using the GNU Programming Environment

The default programming environment is `PrgEnv-cray`, which does not support C++ 11. Switch to the GNU compiler.

```bash
module swap PrgEnv-cray PrgEnv-gnu
```

#### 64-bit Indices with Cray PETSc

Blue Waters provides a `cray-petsc-64` module. However, at the time of writing, it only works with a specific version of the `cray-mpich` library. Use the following to get the module to work:

```
module load cray-petsc-64
module swap cray-mpich cray-mpich/7.2.5
```

### MPI Runner

Blue Waters uses `aprun` instead of `mpirun`.


## Comet

### PETSc Modules Don't Work

The PETSc modules available on Comet were compiled with the Intel 13 compiler. TalyFEM needs to be compiled with 15 or higher. Unfortunately, the Intel 15 runtime is not backwards compatible with the Intel 13 runtime, so you cannot use TalyFEM with the system PETSc modules. To get around this, you'll need to compile PETSc yourself.

### MPI Runner

Comet uses `ibrun` instead of `mpirun`.


## Miscellaneous

### HDF5 Lock Errors on Lustre File Systems

When writing an HDF5 file on a Lustre file system, you might get an error like this:

```
This requires fcntl(2) to be implemented. As of 8/25/2011 it is not. Generic MPICH Message: File locking failed in ADIOI_Set_lock(fd D,cmd F_SETLKW/7,type F_RDLCK/0,whence 0) with return value FFFFFFFF and errno 26.
- If the file system is NFS, you need to use NFS version 3, ensure that the lockd daemon is running on all the machines, and mount the directory with the 'noac' option (no attribute caching).
- If the file system is LUSTRE, ensure that the directory is mounted with the 'flock' option.
ADIOI_Set_lock:: Function not implemented
ADIOI_Set_lock:offset 0, length 1
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
```

My guess for why this happens is that the Intel version of MPI-IO is configured to support file locking, but the Lustre partition is mounted with file locking disabled.

My solution was to use an MPI library with file locking disabled. Adding `--download-mpich` when configuring PETSc seems to fix it. This forces PETSc to download and compile its own version of MPI, which it presumably builds with locking disabled.

### Compiling PETSc on Clusters with a Thread Limit

Some clusters do not allow processes to freely spawn threads (to prevent people from running resource-intensive jobs on the head node, which slows down the system for other users). By default, the PETSc makefile will try to spawn extra processes to speed up build times. If you run into errors related to build threads mysteriously dying, you can add `--useThreads=0` to your configure arguments to disable this.
