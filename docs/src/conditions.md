# Boundary and Initial Conditions

```@index
Pages = ["conditions.md"]
```

```@docs
addBoundaryID(bid::Int, trueOnBdry)
boundary(var, bid, bc_type, bc_exp=0)
referencePoint(var, pos, val)
initial(var, ics)
evalInitialConditions()
```
