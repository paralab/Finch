#=
# 2D NS eq. Dirichlet bc, CG
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("NS");

# Optionally generate a log
useLog("NSlog")

# Set up the configuration (order doesn't matter)
domain(2)                   # dimension
solverType(CG)              # Use CG solver (default)
functionSpace(order=1)      # basis polynomial order
nodeType(LOBATTO)           # GLL elemental node arrangement (default)
timeStepper(EULER_IMPLICIT) # time stepper
dt = 0.05;
nsteps = 50;
setSteps(dt, nsteps)        # manually set stepper dt and Nsteps, overriding defaults and interval

# Specify the problem
mesh(QUADMESH, elsperdim=32, bids=4)

u = variable("u")       #Note: many of these variables are needed for
v = variable("v")       #      linearizing the non-linear problem.
uold = variable("uold") #      Eventually, they will be automatically
vold = variable("vold") #      generated and this will become simpler.
p = variable("p")       #
du = variable("du")     #
dv = variable("dv")     #
dp = variable("dp")     #

testSymbol("w")

# Initial conditions
initial(u, "y > 0.99 ? 1 : 0")
initial(uold, "y > 0.99 ? 1 : 0")
initial(du, "0")
initial(v, "0")
initial(vold, "0")
initial(dv, "0")
initial(p, "0")
initial(dp, "0")

# Boundary conditions
boundary(du, 1, DIRICHLET, 0)
boundary(du, 2, DIRICHLET, 0)
boundary(du, 3, DIRICHLET, 0)
boundary(du, 4, DIRICHLET, 0)

boundary(dv, 1, DIRICHLET, 0)
boundary(dv, 2, DIRICHLET, 0)
boundary(dv, 3, DIRICHLET, 0)
boundary(dv, 4, DIRICHLET, 0)

boundary(dp, 1, NO_BC)
boundary(dp, 2, NO_BC)
boundary(dp, 3, NO_BC)
boundary(dp, 4, NO_BC)
# A single Dirichlet reference point must be used to pin values for dp
referencePoint(dp, [0,0], 0) 

# Write the weak form
coefficient("mu", 0.01)
coefficient("dtc", dt)
coefficient("h", 1.0 / 32)
coefficient("coe1", 4.0)
coefficient("coe2", 6.0)

parameter("tauM", "1.0 ./ (coe1 ./ dtc ./ dtc+ (u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")
parameter("tauC", "0.1*h*h* ((u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")

weakForm([du, dv, dp], ["w*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2) + deriv(u,2)*dv)) - deriv(w,1)*dp + mu*dot(grad(w), grad(du)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1)) - (w*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2))) - deriv(w,1)*p + mu*dot(grad(w), grad(u)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1)))", 
                        "w*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2) + deriv(v,1)*du)) - deriv(w,2)*dp + mu*dot(grad(w), grad(dv)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2)) - (w*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2))) - deriv(w,2)*p + mu*dot(grad(w), grad(v)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2)))", 
	                     "w*(deriv(du,1)+deriv(dv,2)) + tauC*(deriv(w,1)*( du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1) ) + deriv(w,2)*( dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2) )) - (w*(deriv(u,1)+deriv(v,2)) + tauC*(deriv(w,1)*( (u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1) ) + deriv(w,2)*( (v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2) )))"])

#exportCode("NScode")
#importCode("NScode")

solve([du, dv, dp], [u, v, p, uold, vold], nonlinear=true);

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st = :surface, reuse=false))
#display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], v.values[:], st = :surface, reuse=false))

# Dump things to the log if desired
log_dump_config();
log_dump_prob();

finalize_finch() # Finish writing and close any files
