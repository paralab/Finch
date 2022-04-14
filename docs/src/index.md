# Documentation for Finch

Finch is a domain specific language and code generation package for 
partial differential equations. It is intended for people with a good 
understanding of the mathematics of the finite element or finite volume 
methods for solving PDEs. To make full use of Finch some knowledge of 
parallel computing is needed, such as how to launch a multi-process 
instance of Julia using MPI, or launching Julia with multiple threads. 
See [MPI.jl](https://juliaparallel.org/MPI.jl/stable/) and 
[Multi-threading](https://docs.julialang.org/en/v1/base/multi-threading/) 
for refernence.

```@contents
Pages = ["start.md", "configuration.md", "mesh.md", "entities.md", "conditions.md", "equation.md", "solution.md", "reorder.md", "misc.md"]
```
