#=
2D explicit BTE. This assumes MPI is being used
for band parallelization AND CUDA is used.
This version accounts for polarizations.
This partitions the computation amongst both the 
frequency bands and polarizations.
example: 40 bands -> 40 L bands + 15 T bands = 55 possible partitions
=#

### If the Finch package has already been added, use this line #########
#using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

initFinch("FVbte2dgpu");

useLog("FVbte2dgpulog-test", level=3)

# constants and various functions are in another file
includeFile("bte-parameters.jl")
# A set of callback functions for the boundary condition
# includeFile("bte-boundary.jl")
includeFile("bte-boundary-noglobal.jl") # This version does not use global variables
includeFile("bte-post-step.jl")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_EXPLICIT)
dt = 1e-12;
nsteps = 2; # 2 steps for the warmup pass
setSteps(dt, nsteps);

useCUDA(); # This will attempt to use CUDA

# direction and band numbers
ndirs = 20;
frequency_bands = 40;
# ndirs = 4;
# frequency_bands = 10;
(t_bands, l_bands) = get_band_distribution(frequency_bands);
total_bands = t_bands + l_bands;

println("ndirs $(ndirs) total_bands $(total_bands)")

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
mesh(QUADMESH, # quad elements
    elsperdim=[120, 120], # elements in each direction: 20 x 5 uniform grid
    interval=[0, 525e-6, 0, 525e-6],  # interval in each direction
    bids=4, # 4 boundary IDs for this mesh correspond to left, right, bottom, top
    partitions=cell_partitions) # If there are enough procs to also do cell partitioning

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

# Coefficients and related numbers
(dir_x, dir_y) = get_directions_2d(ndirs)
(center_freq, delta_freq, polarizations) = get_band_frequencies(frequency_bands, low_band, high_band);
polarizations_float = [if i == "T" i=0.0 else i=1.0 end for i in polarizations]
group_v = get_group_speeds(center_freq, polarizations);

# These are set as coefficients because they have known values.
Sx = coefficient("Sx", dir_x, type=VAR_ARRAY) # direction x component
Sy = coefficient("Sy", dir_y, type=VAR_ARRAY) # direction y component
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed
cen_freq_coeff = coefficient("center_freq", center_freq, type=VAR_ARRAY)
delta_freq_coeff = coefficient("delta_freq", delta_freq, type=VAR_ARRAY)
polar_coeff = coefficient("polarizations", polarizations_float, type=VAR_ARRAY)
g20xi_coeff = coefficient("g20xi", g20xi, type=VAR_ARRAY)
g20wi_coeff = coefficient("g20wi", g20wi, type=VAR_ARRAY)

# This extra storage is allocated for parallel temperature update
num_cells = length(temperature.values);
uold = zeros(num_cells); # These strange names are taken from the fortran code
gnb = zeros(num_cells);
unew = zeros(num_cells);
uchange = zeros(num_cells);
gna = zeros(num_cells);
uprime = zeros(num_cells);
converged = fill(false, num_cells);

# All isothermal, top central hot spot
boundary(I, 1, FLUX, "isothermal_bdry_noglobal(I, vg, Sx, Sy, center_freq, polarizations, delta_freq, g20xi, g20wi, band, direction, normal, 300, $(ndirs))") # left (ID=1)
boundary(I, 2, FLUX, "isothermal_bdry_noglobal(I, vg, Sx, Sy, center_freq, polarizations, delta_freq, g20xi, g20wi, band, direction, normal, 300, $(ndirs))") # right (ID=2)
boundary(I, 3, FLUX, "isothermal_bdry_noglobal(I, vg, Sx, Sy, center_freq, polarizations, delta_freq, g20xi, g20wi, band, direction, normal, 300, $(ndirs))") # bottom (ID=3)
boundary(I, 4, FLUX, "isothermal_bdry_noglobal(I, vg, Sx, Sy, center_freq, polarizations, delta_freq, g20xi, g20wi, band, direction, normal, (300 + 50*exp(-(x-262e-6)*(x-262e-6)/(5e-9))), $(ndirs))") # top (ID=4)

init_temp = 300; # The initial equilibrium temperature everywhere
initial(I, [equilibrium_intensity(center_freq[b], delta_freq[b], init_temp, polarizations[b]) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq[b], init_temp, polarizations[b]) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], init_temp, polarizations[b]) for b=1:nbands])
initial(temperature, init_temp);
initial(temperatureLast, init_temp);

# assemblyLoops([band, "elements", direction])
assemblyLoops(["elements", band, direction])

# Get integrated intensity from initial equilibrium.
# To get initial values here we have to manually initialize.
evalInitialConditions();
get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

# After each time step the temperature, equilibrium I, and time scales are updated
cell = Indexer("cell", [1, length(temperature.values)])
postStepFunction("update_temperature(temperature, temperatureLast, I, beta, center_freq, delta_freq, polarizations, G_last, G_next, vg, Io, g20xi, g20wi, $ndirs, $nbands, dt, cell)", [cell])

# BTE:
# Dt(int(Iij dx)) = int((Io-Iij)*beta dx) ) - vg * int(Iij * Si.n ds)
conservationForm(I, "(Io[band] - I[direction,band]) * beta[band] + surface(vg[band] * upwind([Sx[direction];Sy[direction]] , I[direction,band]))")

exportCode("bte2dgpucode-new") # uncomment to export generated code to a file
# importCode("bte2dgpucodein") # uncomment to import code from a file

# The fortran version does something odd with volumes
# We need to manually change it here, or adjust the equation
Finch.finch_state.fv_geo_factors.volume .*= 4 * pi / ndirs;

solve(I)

##################################################
# End warmup
# Begin timed run
##################################################

if Finch.finch_state.config.proc_rank == 0
    println("Warm-up result:")
    show(Finch.finch_state.timer_output)
end
Finch.reset_timer!(Finch.finch_state.timer_output)

nsteps = 100;
setSteps(dt, nsteps);

evalInitialConditions();
get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

solve(I)

# outputValues(temperature, "bte2dTemp", format="vtk");
outputValues(temperature, "bte2dTemp-new", format="csv");

finalizeFinch()

### uncomment below to plot

xy = Finch.finch_state.fv_info.cellCenters
using Plots
pyplot();
display(plot(xy[1,:], xy[2,:], temperature.values[:], st=:surface))
