# Other Interface Functions

These miscellaneous parts are part of the interface, but don't fit in the other categories.

```@index
Pages = ["misc.md"]
```

```@docs
initFinch(name="unnamedProject")
customOperator(name, handle)
customOperatorFile(filename)
VariableTransform
variableTransform(var1, var2, func)
transformVariable(xform::VariableTransform)
preStepFunction(fun)
postStepFunction(fun)
callbackFunction(fun; name="", args=[], body="")
assemblyLoops(indices, parallel_type=[])
exportCode(filename)
importCode(filename)
finalizeFinch()
cachesim(use)
cachesimSolve(var)
```
