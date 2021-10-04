#=
# Linear elasticity
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

init_finch("elasticitynew");

# Optionally generate a log
useLog("elasticitymatlablog", level=3)
generateFor("target_matlab_cg.jl")

n = [30,6,6]; # number of elements in x,y,z
interval = [0,1,0,0.2,0,0.2]; # domain bounds

# Set up the configuration
domain(3)
functionSpace(order=2)

# Specify the problem
mesh(HEXMESH, elsperdim=n, bids=4, interval=interval)

u = variable("u", VECTOR)          # make a vector variable with symbol u
testSymbol("v", VECTOR)            # sets the symbol for a test function

boundary(u, 1, DIRICHLET, [0,0,0]) # x=0
boundary(u, 2, NEUMANN, [0,0,0])   # elsewhere
boundary(u, 3, NEUMANN, [0,0,0])
boundary(u, 4, NEUMANN, [0,0,0])

# Write the weak form
# coefficient("mu", "x>0.5 ? 0.2 : 10") # discontinuous mu
coefficient("mu", 1) # constant mu
coefficient("lambda", 1.25)
coefficient("f", [0,0,-10], VECTOR)
weakForm(u, "inner( (lambda * div(u) .* [1 0 0; 0 1 0; 0 0 1] + mu .* (grad(u) + transpose(grad(u)))), grad(v)) - dot(f,v)")

exportCode("elasticitycode")

solve(u);

finalize_finch() # Finish writing and close any files
