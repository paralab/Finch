#=
3D BTE with GPU
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
dt = 1e-12;
nsteps = 1;
setSteps(dt, nsteps);

useCUDA(); # This will attempt to use CUDA

# direction and band numbers
nphi = 8; # number of phi angles
ntheta = 8; # number of theta angles
ndirs = nphi*ntheta; # total number of directions
frequency_bands = 40; # Set to 40 for larger scale computation
(t_bands, l_bands) = get_band_distribution(frequency_bands);
total_bands = t_bands + l_bands;

# This manual partitioning is needed because this problem is complicated.
np = Finch.finch_state.config.num_procs;
rank = Finch.finch_state.config.proc_rank;
cell_partitions = 1; # only band parallel
num_band_partitions = np;
# Divide bands among processes for band based parallel
if np > total_bands 
    println("Too many processes. Max is "*string(total_bands));
    exit(1);
else
    nbands = Int(floor(total_bands / np));
    low_band = rank*nbands + 1;
    high_band = (rank+1)*nbands;
    extras = mod(total_bands, np);
    if extras > rank
        nbands += 1;
        low_band += rank;
        high_band += rank+1;
    else
        low_band += extras;
        high_band += extras;
    end
end
MPI = Finch.MPI;

# A simple mesh is internally generated for convenience
# This matches the mesh in model_setup_BTE.in
# mesh(HEXMESH, # hex elements
#     elsperdim=[10,10,10], # elements in each direction in uniform grid
#     interval=[0, 10e-6, 0, 10e-6, 0, 10e-6],  # interval in each direction
#     bids=4) # 4 boundary IDs for this mesh correspond to x=0,x=1,y,z sets of faces
    
mesh("cubecorner15k.msh")

addBoundaryID(2, (x,y,z) -> (z < (0.000525 - 1e-15))); # everywhere but the top
addBoundaryID(3, (x,y,z) -> (y < 1e-16 || x < 1e-16)); # the x=0 and y=0 faces

# Indices, Variables, etc.
direction = index("direction", range = [1,ndirs])
band = index("band", range = [1,nbands])

# These are all set as variables because they are unknown, but only I is solved for in the PDE.
I = variable("I", type=VAR_ARRAY, location=CELL, index = [direction, band]) # Intensity
Io = variable("Io", type=VAR_ARRAY, location=CELL, index = [band]) # Equilibrium intensity for each band
beta = variable("beta", type=VAR_ARRAY, location=CELL, index = [band]) # Relaxation time scale
temperature = variable("temperature", location=CELL) # temperature of each cell
temperatureLast = variable("temperatureLast", location=CELL) # temperature from last time step
G_last = variable("G_last", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity from last step
G_next = variable("G_next", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity for current step

bdrything = variable("bdrything", location=CELL)

# Coefficients and related numbers
(dir_x, dir_y, dir_z, omega) = get_directions_3d(nphi, ntheta);
# reflect = get_3d_cube_reflections(dir_x, dir_y, dir_z);
(center_freq, delta_freq, polarizations) = get_band_frequencies(frequency_bands, low_band, high_band);
group_v = get_group_speeds(center_freq, polarizations);

# These are set as coefficients because they have known values.
Sx = coefficient("Sx", dir_x, type=VAR_ARRAY) # direction x component
Sy = coefficient("Sy", dir_y, type=VAR_ARRAY) # direction y component
Sz = coefficient("Sz", dir_z, type=VAR_ARRAY) # direction z component
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed

num_cells = length(temperature.values);
uold = zeros(num_cells); # These strange names are taken from the fortran code
gnb = zeros(num_cells);
unew = zeros(num_cells);
uchange = zeros(num_cells);
gna = zeros(num_cells);
uprime = zeros(num_cells);
converged = fill(false, num_cells);

assemblyLoops(["elements", band, direction])

# boundary(I, 1, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 350)") # (300 + 50*sin(pi*y/(5e-6))*sin(pi*z/(5e-6)))") # x=0
# boundary(I, 2, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 300)") # x=1
# boundary(I, 3, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 300)") # y faces
# boundary(I, 4, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 300)") # z faces

boundary(I, 1, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, (300 + 50*exp(-(x*x+y*y)/(0.000000005))))")
boundary(I, 2, FLUX, "isothermal_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal, 300)")
boundary(I, 3, FLUX, "symmetric_bdry_3d(I, vg, Sx, Sy, Sz, band, direction, normal)")

init_temp = 300; # The initial equilibrium temperature everywhere
initial(I, [equilibrium_intensity(center_freq[b], delta_freq[b], init_temp, polarizations[b]) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq[b], init_temp, polarizations[b]) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], init_temp, polarizations[b]) for b=1:nbands])
initial(temperature, init_temp);
initial(temperatureLast, init_temp);

# To get initial values here before calling solve, manually initialize.
evalInitialConditions();
get_integrated_intensity_3d!(G_last.values, I.values, ndirs, nbands, omega);

# After each time step the temperature, equilibrium I, and time scales are updated
function post_step()
    update_temperature(temperature.values, temperatureLast.values, I.values, beta.values, 
                        center_freq, delta_freq, converged, polarizations, band_parallel=true, omega=omega);
end
postStepFunction(post_step);

# BTE:
# Dt(int(Iij dx)) = int((Io-I)/beta dx) ) + vg * int(I * S.n ds)
# Input conservation form representing: (Io-I)/beta + surface(vg * I * S.n)
conservationForm(I, "(Io[band] - I[direction,band]) * beta[band] + surface(vg[band] * upwind([Sx[direction];Sy[direction];Sz[direction]] , I[direction,band]))")

exportCode("bte3dgpucode") # uncomment to export generated code to a file
# importCode("bte3dgpucodein") # uncomment to import code from a file

solve(I)

outputValues(temperature, "bte3dgpuTemp", format="vtk");

finalizeFinch()
