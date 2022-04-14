# Entities

Entities include variables, coefficients, indices, and other objects that
appear in the equations.

```@index
Pages = ["entities.md"]
```

```@docs
Variable
Coefficient
Indexer
variable(name; type=SCALAR, location=NODAL, method=CG, index=nothing)
coefficient(name, val; type=SCALAR, location=NODAL, element_array=false)
parameter(name, val; type=SCALAR)
index(name; range=[1])
testSymbol(symbol; type=SCALAR)
```
