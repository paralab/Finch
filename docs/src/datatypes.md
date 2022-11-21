# Data types

The data type used for many of the integer and floating point arrays 
can be set when initializing Finch. The floating point type must be a 
subtype of AbstractFloat. Most of the arrays used in computation, such 
as variable value arrays and the global linear system, will hold this 
type of data.

To set a type other than the defalt `Float64`, pass the desired type
to `initFinch("projectName", myFloatType)`.

Note that some parts of the computation may involve default types like 
Float64 and Int, so mixed arithmetic operations and conversions should be 
defined. In particular, Finch's use of MPI is not yet configured to 
communicate custom data types for which `isbitstype(T)` is `false`, 
so unexpected conversions may occur when using MPI.

```@index
Pages = ["datatypes.md"]
```