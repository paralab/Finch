#=
# 2D Poisson, Dirichlet bc
# Uses Matlab imported as a custom gen target
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("poisson2dcustomnew");

# Try making an optional log
useLog("poisson2dcustomnewlog")

# Generate for the target in customtarget.jl
generateFor("target_matlab_cg.jl")

# Set up the configuration (order doesn't matter)
domain(2)
functionSpace(order=2)

# Specify the problem
mesh(QUADMESH, elsperdim=30)

u = variable("u")

testSymbol("v")

boundary(u, 1, DIRICHLET, "sin(3*pi*x)")

# Write the weak form 
coefficient("f", "-2*pi*pi*sin(pi*x)*sin(pi*y)")
weakForm(u, "-dot(grad(u),grad(v)) - f*v")

solve(u);

finalize_finch()
