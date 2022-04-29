# Heat

```@raw html
<img src="../assets/images/heat.png" alt="heat" width="400">
```

The script file: [example-heat2d.jl](https://github.com/paralab/Finch/blob/master/Finch/examples/example-heat2d.jl)

A 2D heat equation demonstrates support for time dependent problems.

Begin by importing and using the Finch module. Then initialize. The name here is only used when generating code files.
```
using Finch
init_finch("heat2d");
```
Then set up the configuration. This example simply sets dimensionality of the domain and polynomial order of the basis function space.
```
domain(2)                  	# dimension
functionSpace(order=4) 		# polynomial order
```
Use the built-in simple mesh generator to make the mesh and set up all node mappings.
```
mesh(QUADMESH, elsperdim=10)# has 10*10 uniform, square elements
```
Define the variable, test function, and coefficient symbols.
```
u = variable("u")           # make a scalar variable u
testSymbol("v")             # sets the symbol for a test function

coefficient("f", "0.5*sin(6*pi*x)*sin(6*pi*y)")
```
Set up the time stepper and initial conditions. This example uses a low-storage RK4. Other explicit or implicit methods are available.
```
timeStepper(LSRK4)  		# Low-storage RK4
timeInterval(1) 			# The end time
initial(u, "abs(x-0.5)+abs(y-0.5) < 0.2 ? 1 : 0") # initial condition
```
Convert the PDE
```@raw html
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=\frac{d}{dt}u%2BD\Delta%20u=f"> </div>
```
into the weak form
```@raw html
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=\frac{d}{dt}(u,v)%2BD(\nabla%20u,\nabla%20v)=(f,v)"> </div>
```

The boundary conditions are specified.
```
boundary(u, 1, DIRICHLET, 0)
```
Then write the weak form expression in the residual form. Finally, solve for u.
```
weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")
solve(u);
```
End things with `finalize_finch()` to finish up any generated files and the log.