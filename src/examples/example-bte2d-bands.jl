#=
2D explicit BTE. This assumes MPI is being used
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("FVbte2dband");

useLog("FVbte2dbandlog", level=3)

# constants and various functions are in another file
include("bte-parameters.jl")
# A set of callback functions for the boundary condition
include("bte-boundary.jl")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_EXPLICIT)
dt = 2e-12;
nsteps = 2; # only 2 warm-up steps. More steps for final run.
setSteps(dt, nsteps);

# direction and band numbers
ndirs = 16;
total_bands = 4;
np = Finch.finch_state.config.num_procs;
rank = Finch.finch_state.config.proc_rank;
cell_partitions = 1; # This will only be increased if at least 2*total_bands procs
# Divide bands among processes for primarily band based parallel
if np >= total_bands # If there are extra processes, use them for cell partitioning
    cell_partitions = Int(floor(np / total_bands));

    nbands = 1;
    low_band = Int(floor(rank/cell_partitions)) + 1;
    high_band = low_band;

else
    nbands = Int(floor(total_bands / np));
    low_band = rank*nbands + 1;
    high_band = (rank+1)*nbands;
    if mod(total_bands, np) > rank
        nbands += 1;
        low_band += rank;
        high_band += rank+1;
    end
end
MPI = Finch.MPI;
if cell_partitions > 1
    partition_num = Int(floor((rank+0.5) * cell_partitions / Finch.finch_state.config.num_procs));
    cell_comm = MPI.Comm_split(MPI.COMM_WORLD, partition_num, rank); # communicator for this mesh partition
    # cell_rank = MPI.Comm_rank(cell_comm);
else
    cell_comm = MPI.Comm_dup(MPI.COMM_WORLD);
end
num_band_partitions = MPI.Comm_size(cell_comm);
if num_band_partitions > total_bands
    println("Error: number of band partitions = "*string(num_band_partitions)*", but there are only "*string(total_bands)*" bands.");
    exit(1);
end

# A simple mesh is internally generated for convenience
# This matches the mesh in model_setup_BTE.in
mesh(QUADMESH, # quad elements
    elsperdim=[20,20], # elements in each direction: 20 x 5 uniform grid
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
(dir_x, dir_y, reflect) = get_directions_2d(ndirs)
(center_freq, delta_freq) = get_band_frequencies(total_bands, low_band, high_band);
group_v = get_group_speeds(center_freq);

# These are set as coefficients because they have known values.
Sx = coefficient("Sx", dir_x, type=VAR_ARRAY) # direction x component
Sy = coefficient("Sy", dir_y, type=VAR_ARRAY) # direction y component
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed

# This extra storage is allocated for parallel temperature update
num_cells = length(temperature.values);
uold = zeros(num_cells*num_band_partitions); # These strange names are taken from the fortran code
gnb = zeros(num_cells*num_band_partitions);
unew = zeros(num_cells*num_band_partitions);
gna = zeros(num_cells*num_band_partitions);
uprime = zeros(num_cells*num_band_partitions);

uold_buffer = MPI.UBuffer(uold, num_cells, num_band_partitions, MPI.Datatype(Float64));
gnb_buffer = MPI.UBuffer(gnb, num_cells, num_band_partitions, MPI.Datatype(Float64));
unew_buffer = MPI.UBuffer(unew, num_cells, num_band_partitions, MPI.Datatype(Float64));
gna_buffer = MPI.UBuffer(gna, num_cells, num_band_partitions, MPI.Datatype(Float64));
uprime_buffer = MPI.UBuffer(uprime, num_cells, num_band_partitions, MPI.Datatype(Float64));

boundary(I, 1, FLUX, "symmetric_bdry(I, vg, Sx, Sy, band, direction, normal)") # left (ID=1)
# boundary(I, 2, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, 300)") # right (ID=2)
boundary(I, 2, FLUX, "symmetric_bdry(I, vg, Sx, Sy, band, direction, normal)") # left (ID=1)
boundary(I, 3, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, 300)") # bottom (ID=3)
# boundary(I, 4, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, (300 + 10*exp(-x*x/(32e-12))))") # top (ID=4)
boundary(I, 4, FLUX, "isothermal_bdry(I, vg, Sx, Sy, band, direction, normal, 305)") # top (ID=4)

init_temp = 300; # The initial equilibrium temperature everywhere
initial(I, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq, init_temp) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], init_temp) for b=1:nbands])
initial(temperature, init_temp);
initial(temperatureLast, init_temp);

# assemblyLoops([band, "elements", direction])
assemblyLoops(["elements", band, direction])

# Create an array to hold the values of I from the last step.
# To get initial values here we have to manually initialize.
evalInitialConditions();
get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

# After each time step the temperature, equilibrium I, and time scales are updated
function post_step()
    update_temperature(temperature.values, temperatureLast.values, I.values, beta.values, center_freq, delta_freq);
end
postStepFunction(post_step);

# BTE:
# Dt(int(Iij dx)) = int((Io-Iij)*beta dx) ) - vg * int(Iij * Si.n ds)
conservationForm(I, "(Io[band] - I[direction,band]) / beta[band] + surface(vg[band] * upwind([Sx[direction];Sy[direction]] , I[direction,band]))")

exportCode("bte2dcode") # uncomment to export generated code to a file
# importCode("bte2dcode") # uncomment to import code from a file

solve(I)

if Finch.finch_state.config.proc_rank == 0
    println("Warm-up result:")
    show(Finch.finch_state.timer_output)
end
Finch.reset_timer!(Finch.finch_state.timer_output)

##################################################
# End warmup, start real
##################################################
# nsteps = 10000;
# setSteps(dt, nsteps);

# evalInitialConditions();
# get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

# solve(I)

outputValues(temperature, "bte2dTemp", format="vtk");

finalizeFinch()
