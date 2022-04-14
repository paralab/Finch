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
matrixFree(shallwe=true; maxiters=100, tol=1e-6)
```
