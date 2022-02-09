#=
# 1D mixed FE/FV example
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("fefv1d");

useLog("fefv1dlog")

# Configuration setup
domain(1)
solverType([CG,FV])
timeStepper(EULER_EXPLICIT)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n)

# Variables and BCs
u = variable("u", method=CG)
testSymbol("v")
boundary(u, 1, DIRICHLET, 0)

w = variable("w", SCALAR, CELL, method=FV)
boundary(w, 1, DIRICHLET, 0)

# Time interval and initial condition
T = 1;
timeInterval(T)
initial(u, "x<0.5 ? sin(2*pi*x)^4 : 0")
initial(w, "x>0.5 ? sin(2*pi*x)^4 : 0")

coefficient("D", 0.001) # Diffusion rate
coefficient("f", "-0.5*sin(pi*x)^2")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F ds)
flux(w, "-D*dot(grad(w),normal())") # in 1D -> -D * dx(u) * normal
source(w, "u")

# The weak form expression
weakForm(u, "Dt(u*v) + D * dot(grad(u),grad(v)) - w*v")

solve([u,w])

finalize_finch()

##### Uncomment below to plot
xc = Finch.fv_info.cellCenters[:];
xn = Finch.grid_data.allnodes[1,:];
nnodes = length(xn); # The initial condition
u0 = zeros(nnodes);
for i=1:nnodes
    u0[i] = sin(2*pi*xn[i])^4
end

using Plots
pyplot();
display(plot([xn xn], [u0 u.values[:]], markershape=:circle, label=["initial" "u, t="*string(T)]))
plot!(xc, w.values[:], markershape=:square, label="w, t="*string(T))
