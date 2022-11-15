# Configuration

These functions are used at the beginning of your script to set up 
the configuration.

```@index
Pages = ["configuration.md"]
```

```@docs
generateFor(lang; filename=project_name, header="", params=nothing)
useLog(name=project_name; dir=output_dir, level=2)
domain(dims; shape=SQUARE, grid=UNIFORM_GRID)
solverType(method, backend=DEFAULT_SOLVER)
functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
nodeType(type)
timeStepper(type; cfl=0)
timeInterval(T)
setSteps(dt, steps)
linAlgOptions(;matrixFree::Bool=false, iterative::Bool=false, method::String="GMRES", pc::String="ILU", maxiter::Int=0, abstol=0, reltol=1e-8, gmresRestart::Int=0, verbose::Bool=false)
```
