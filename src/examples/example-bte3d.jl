#=
3D BTE
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("FVbte3d");
useLog("FVbte3dlog", level=3)

# constants and various functions are in another file
include("bte-parameters.jl")
# A set of callback functions for the boundary condition
include("bte-boundary.jl")

# Configuration setup
domain(3)
solverType(FV)
timeStepper(EULER_EXPLICIT)
# Specify time steps if desired
dt = 2.5e-12;
nsteps = 10;
setSteps(dt, nsteps);

# direction and band numbers
nphi = 8; # number of phi angles
ntheta = 8; # number of theta angles
ndirs = nphi*ntheta; # total number of directions
nbands = 4; # Set to 40 for larger scale computation

# A simple mesh is internally generated for convenience
# This matches the mesh in model_setup_BTE.in
mesh(HEXMESH, # hex elements
    elsperdim=[10,5,5], # elements in each direction in uniform grid
    interval=[0, 4e-6, 0, 2e-6, 0, 2e-6, 0, 2e-6],  # interval in each direction
    bids=4) # 4 boundary IDs for this mesh correspond to x=0,x=1,y,z sets of faces

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
(dir_x, dir_y, dir_z, omega) = get_directions_3d(nphi, ntheta);
reflect = get_3d_cube_reflections(dir_x, dir_y, dir_z);
(center_freq, delta_freq) = get_band_frequencies(nbands);
group_v = get_group_speeds(center_freq);

# These are set as coefficients because they have known values.
Sx = coefficient("Sx", dir_x, type=VAR_ARRAY) # direction x component
Sy = coefficient("Sy", dir_y, type=VAR_ARRAY) # direction y component
Sz = coefficient("Sz", dir_z, type=VAR_ARRAY) # direction z component
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed

boundary(I, 1, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 305)") # x=0
boundary(I, 2, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 300)") # x=1
boundary(I, 3, FLUX, "symmetric_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal)") # y faces
boundary(I, 4, FLUX, "symmetric_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal)") # z faces

init_temp = 300; # The initial equilibrium temperature everywhere
initial(I, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], init_temp) for b=1:nbands])
initial(temperature, init_temp);

# To get initial values here before calling solve, manually initialize.
evalInitialConditions();
get_integrated_intensity_3d!(G_last.values, I.values, ndirs, nbands, omega);

# After each time step the temperature, equilibrium I, and time scales are updated
function post_step()
    update_temperature(temperature.values, I.values, center_freq, delta_freq, threed=true, omega=omega);
end
postStepFunction(post_step);

assemblyLoops(["elements", band, direction])

# BTE:
# Dt(int(Iij dx)) = int((Io-I)/beta dx) ) + vg * int(I * S.n ds)
# Input conservation form representing: (Io-I)/beta + surface(vg * I * S.n)
conservationForm(I, "(Io[band] - I[direction,band]) / beta[band] + surface(vg[band] * upwind([Sx[direction];Sy[direction];Sz[direction]] , I[direction,band]))")

exportCode("bte3dcode") # uncomment to export generated code to a file
# importCode("bte3dcode") # uncomment to import code from a file

solve(I)

output_values(temperature, "bte3dTemp", format="vtk");

finalizeFinch()
