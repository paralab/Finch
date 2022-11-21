#=
2D BTE
Can use EULER_EXPLICIT or EULER_IMPLICIT steppers.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("FVbte2d");
useLog("FVbte2dlog", level=3)

# constants and various functions are in another file
include("bte-parameters.jl")
# A set of callback functions for the boundary condition
include("bte-boundary.jl")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_EXPLICIT)
# Specify time steps if desired
dt = 2.5e-11;
nsteps = 100;
setSteps(dt, nsteps);

# direction and band numbers
ndirs = 16;
nbands = 4; # Set to 40 for larger scale computation

# A simple mesh is internally generated for convenience
# This matches the mesh in model_setup_BTE.in
mesh(QUADMESH, # quad elements
    elsperdim=[20,5], # elements in each direction: 20 x 5 uniform grid
    interval=[0, 4e-6, 0, 1e-6],  # interval in each direction: a very small rectangle
    bids=3) # 3 boundary IDs for this mesh correspond to left, right, top/bottom

# Indices, Variables, etc.
direction = index("direction", range = [1,ndirs])
band = index("band", range = [1,nbands])

# These are all set as variables because they are unknown, but only I is solved for in the PDE.
I = variable("I", type=VAR_ARRAY, location=CELL, index = [direction, band]) # Intensity
Io = variable("Io", type=VAR_ARRAY, location=CELL, index = [band]) # Equilibrium intensity for each band
beta = variable("beta", type=VAR_ARRAY, location=CELL, index = [band]) # Relaxation time scale
temperature = variable("temperature", location=CELL) # temperature of each cell
G_last = variable("G_last", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity from last step
G_next = variable("G_next", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity for current step

# Coefficients and related numbers
(dir_x, dir_y, reflect) = get_directions_2d(ndirs)
(center_freq, delta_freq) = get_band_frequencies(nbands);
group_v = get_group_speeds(center_freq);

# These are set as coefficients because they have known values.
Sx = coefficient("Sx", dir_x, type=VAR_ARRAY) # direction x component
Sy = coefficient("Sy", dir_y, type=VAR_ARRAY) # direction y component
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed

boundary(I, 1, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, 305)") # left (ID=1)
boundary(I, 2, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, 300)") # right (ID=2)
boundary(I, 3, FLUX, "symmetric_bdry(I, vg, Sx, Sy, band, direction, normal)") # top and bottom (ID=3)

init_temp = 300; # The initial equilibrium temperature everywhere
initial(I, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], init_temp) for b=1:nbands])
initial(temperature, init_temp);

# To get initial values here before calling solve, manually initialize.
evalInitialConditions();
get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

# After each time step the temperature, equilibrium I, and time scales are updated
function post_step()
    update_temperature(temperature.values, I.values, center_freq, delta_freq);
end
postStepFunction(post_step);

assemblyLoops(["elements", band, direction])

# BTE:
# Dt(int(Iij dx)) = int((Io-I)/beta dx) ) + vg * int(I * S.n ds)
# Input conservation form representing: (Io-I)/beta + surface(vg * I * S.n)
conservationForm(I, "(Io[band] - I[direction,band]) / beta[band] + surface(vg[band] * upwind([Sx[direction];Sy[direction]] , I[direction,band]))")

exportCode("bte2dcode") # uncomment to export generated code to a file
# importCode("bte2dcode") # uncomment to import code from a file

solve(I)

outputValues(temperature, "bte2dTemp", format="vtk");

finalizeFinch()

##### Uncomment below to plot ######

# using Plots
# pyplot();

# xy = Finch.fv_info.cellCenters

# # plots the temperature
# display(plot(xy[1,:], xy[2,:], temperature.values[:], st=:surface));

# # plots the intensity for the first band an 16 directions
# p1 = plot(xy[1,:], xy[2,:], I.values[1,:], st=:surface)#, zlims=(0,Inf))
# p2 = plot(xy[1,:], xy[2,:], I.values[2,:], st=:surface)#, zlims=(0,Inf))
# p3 = plot(xy[1,:], xy[2,:], I.values[3,:], st=:surface)#, zlims=(0,Inf))
# p4 = plot(xy[1,:], xy[2,:], I.values[4,:], st=:surface)#, zlims=(0,Inf))
# p5 = plot(xy[1,:], xy[2,:], I.values[5,:], st=:surface)#, zlims=(0,Inf))
# p6 = plot(xy[1,:], xy[2,:], I.values[6,:], st=:surface)#, zlims=(0,Inf))
# p7 = plot(xy[1,:], xy[2,:], I.values[7,:], st=:surface)#, zlims=(0,Inf))
# p8 = plot(xy[1,:], xy[2,:], I.values[8,:], st=:surface)#, zlims=(0,Inf))
# p9 = plot(xy[1,:], xy[2,:], I.values[9,:], st=:surface)#, zlims=(0,Inf))
# p10 = plot(xy[1,:], xy[2,:], I.values[10,:], st=:surface)#, zlims=(0,Inf))
# p11 = plot(xy[1,:], xy[2,:], I.values[11,:], st=:surface)#, zlims=(0,Inf))
# p12 = plot(xy[1,:], xy[2,:], I.values[12,:], st=:surface)#, zlims=(0,Inf))
# p13 = plot(xy[1,:], xy[2,:], I.values[13,:], st=:surface)#, zlims=(0,Inf))
# p14 = plot(xy[1,:], xy[2,:], I.values[14,:], st=:surface)#, zlims=(0,Inf))
# p15 = plot(xy[1,:], xy[2,:], I.values[15,:], st=:surface)#, zlims=(0,Inf))
# p16 = plot(xy[1,:], xy[2,:], I.values[16,:], st=:surface)#, zlims=(0,Inf))
# display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, layout=16))
