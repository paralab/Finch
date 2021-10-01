if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("FVadvection1d");

useLog("FVadvection1dlog")

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n, bids=2)

# Variables and BCs
u = variable("u", SCALAR, CELL)
v = variable("v", SCALAR, CELL)
boundary(u, 1, FLUX, "t<0.2 ? 1 : 0")
boundary(u, 2, NO_BC)
boundary(v, 1, FLUX, "t<0.2 ? 1 : 0")
boundary(v, 2, NO_BC)

# Time interval and initial condition
T = 0.5;
timeInterval(T)
initial(u, "0")
initial(v, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("a", 1) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
flux([u, v], ["upwind(a,u)", "upwind(a,v,0.75)"]) 
# Note that there is no source() for this problem. 

#@exportCode("fvad1dcode") # uncomment to export generated code to a file
#@importCode("fvad1dcode") # uncomment to import code from a file

solve([u,v])

finalize_finch()

##### Uncomment below to compare to exact solution

# # The exact solution with constant velocity v
# a = 1;
# exact = zeros(n);
# x = Finch.fv_info.cellCenters[:]
# for i=1:n
#     xt = x[i] - a*T;
#     if xt < 0
#         exact[i] = xt > -0.2 ? 1 : 0
#     end
# end

# using Plots
# pyplot();
# display(plot([x x x], [exact u.values[:] v.values[:]], markershape=:circle, label=["exact" "upwind" "alpha=0.75"]))
