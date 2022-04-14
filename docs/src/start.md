# Getting started

Add the Finch package using `]add "https://github.com/paralab/Finch"`

Alternatively, you can use a local copy of the source by navigating to 
the directory containing the Finch.jl file and typing 
`include("Finch.jl"); using .Finch`

Typically a Julia script file will be created that will perform these 
tasks in this order.
1. Set up configuration.
2. Create or import mesh.
3. Define variables and other entities.
4. Add boundary and initial conditions.
5. Input the PDE expressions.
6. Solve.
7. Process or output the results.

See the example scripts for a more detailed illustration of the structure.
