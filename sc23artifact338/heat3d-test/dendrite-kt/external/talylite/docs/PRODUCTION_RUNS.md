# Production Run Checklist

### Are you using an optimized build of PETSc?

By default, PETSc is built with debugging turned on. To turn it off, add `--with-debugging=0` to your `./configure` options.

### Are you using an optimized build of TalyFEM?

Configure CMake with `cmake .. -DCMAKE_BUILD_TYPE=Release` or `-DCMAKE_BUILD_TYPE=RelWithDebInfo`.

### Have you enabled architecture-specific instruction sets?

You can add `-xHost` or `-mtune=native` to your compiler options to enable all architecture-specific instruction sets available on the processor you are compiling code on (e.g. SSE, AVX). Be careful on clusters where the head node's processor does not match the compute node's processor.

### Do you need 64-bit PETSc indices?

If your problem will have more than 2,147,483,647 total degrees of freedom in your problem (rows in the matrix), you will need to compile PETSc with 64-bit PETSc indices. To do this, configure PETSc with `--with-64-bit-indices` and recompile TalyFEM.

### Have you tested everything?

Do a test run with a smaller problem first to make sure your configuration is working.

### Do you need reproducibility information?

TalyFEM can generate a file for you that includes information about compile-time and run-time information for later auditing. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for more details. This is less useful if you are running thousands of smaller jobs.
