#=
# 2D NS eq. Dirichlet bc, CG
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("NSu");

@useLog("NSulog")
set_log_level(3)

# Set up the configuration
@domain(2, IRREGULAR, UNSTRUCTURED)
@solver(CG)
@functionSpace(LEGENDRE, 2)
@nodes(LOBATTO)
@stepper(EULER_IMPLICIT)
dt = 0.05;
nsteps = 50;
@setSteps(dt, nsteps)  # manually set stepper dt and Nsteps, overriding defaults and interval

# Specify the problem
@mesh("src/examples/rect.msh") # this rectangle covers [0, 0.1]x[0, 0.3]

# left and right boundaries have bid 2, 3
@addBoundaryID(2, x <= 0)
@addBoundaryID(3, x >= 0.1)

@variable(u)
@variable(v)
@variable(uold)
@variable(vold)
@variable(p)
@variable(du)
@variable(dv)
@variable(dp) 

@testSymbol(w)

#T = 2 # set manually using setSteps above
#@timeInterval(T)
@initial(u, "0")
@initial(uold, "0")
@initial(du, "0")
@initial(v, "(x < 0.001 || x > 0.099) ? -1*sign(x-1) : 0") # up on left, down on right
@initial(vold, "(x < 0.001 || x > 0.099) ? -1*sign(x-1) : 0")
@initial(dv, "0")
@initial(p, "0")
@initial(dp, "0")

@boundary(du, 1, DIRICHLET, 0)
@boundary(du, 2, DIRICHLET, 0)
@boundary(du, 3, DIRICHLET, 0)
#@boundary(du, 4, DIRICHLET, 0)

@boundary(dv, 1, DIRICHLET, 0)
@boundary(dv, 2, DIRICHLET, 0)
@boundary(dv, 3, DIRICHLET, 0)
#@boundary(dv, 4, DIRICHLET, 0)

@boundary(dp, 1, NO_BC)
@boundary(dp, 2, NO_BC)
@boundary(dp, 3, NO_BC)
#@boundary(dp, 4, NO_BC)
@referencePoint(dp, [0.05,0], 0)


# Write the weak form
@coefficient(mu, 0.01)
@coefficient(dtc, dt)
@coefficient(h, 0.01)
@coefficient(coe1, 4.0)
@coefficient(coe2, 6.0)

@parameter(tauM, "1.0 ./ (coe1 ./ dtc ./ dtc+ (u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")
@parameter(tauC, "0.1*h*h* ((u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")

@weakForm([du, dv, dp], ["w*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2) + deriv(u,2)*dv)) - deriv(w,1)*dp + mu*dot(grad(w), grad(du)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1)) - (w*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2))) - deriv(w,1)*p + mu*dot(grad(w), grad(u)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1)))", 
                         "w*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2) + deriv(v,1)*du)) - deriv(w,2)*dp + mu*dot(grad(w), grad(dv)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2)) - (w*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2))) - deriv(w,2)*p + mu*dot(grad(w), grad(v)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2)))", 
	                     "w*(deriv(du,1)+deriv(dv,2)) + tauC*(deriv(w,1)*( du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1) ) + deriv(w,2)*( dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2) )) - (w*(deriv(u,1)+deriv(v,2)) + tauC*(deriv(w,1)*( (u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1) ) + deriv(w,2)*( (v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2) )))"])

@exportCode("NScode")

solve([du, dv, dp], [u, v, p, uold, vold], nonlinear=true);

# output to file
output_values([u,v], "NSunstructured", format="vtk");

# u and v values are written to files, but the node locations are needed
# g = Finch.grid_data.allnodes;
# N = size(g, 2);
# outfile = "grid_xy";
# open(outfile, "w") do f
#     for i=1:N
#         println(f, string(g[1,i]));
#         println(f, string(g[2,i]));
#     end
#     close(f)
# end

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st = :surface, reuse=false))
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], v.values[:], st = :surface, reuse=false))

# check
# log_dump_config(Finch.config);
# log_dump_prob(Finch.prob);

@finalize()
