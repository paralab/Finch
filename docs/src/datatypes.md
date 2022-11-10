# Data types

The data type used for many of the integer and floating point arrays 
can be set. The floating point type is a subtype of AbstractFloat. 
Most of the arrays used in computation, such as variable value arrays 
and the global linear system, will hold this type of data. 

The index data type is a subtype of the abstract Integer 
type. It is used for most index maps such as element connectivity maps.

Note that some parts of the computation may involve default types like 
Float64 and Int, so mixed arithmetic operations and conversions should be 
defined. In particular, Finch's use of MPI is not yet configured to 
communicate custom data types, so unexpected conversions may occur 
when using MPI.

When using these optional commands, be sure to call them near the beginning 
of the script before setting a mesh or creating any variables. That 
way things will be initialized with the desired type.

```@index
Pages = ["datatypes.md"]
```

```@docs
floatDataType(type)
indexDataType(type)
```
