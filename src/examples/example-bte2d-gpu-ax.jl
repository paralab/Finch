#=
2D axisymmetric explicit BTE with transducer(ALSI) BC
Uses MPI to do bands in parallel.
This will also use CUDA.

This version accounts for polarizations as separate bands.
This partitions the computation amongst both the 
frequency bands and polarizations.
example: 40 bands -> 40 L bands + 15 T bands = 55 possible partitions
=#

### If the Finch package has already been added, use this line ###########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("bteAxiBand");
useLog("bteAxiBandlog", level=3)

# constants and various functions are in another file
include("bte-parameters-axi.jl")
# A set of callback functions for the boundary conditions
include("bte-boundary-axi.jl")

# Input parameters ###########################################################
nphi = 4;           # number of phi angles
ntheta = 20;         # number of theta angles
frequency_bands = 40;# number of frequency bands (not including polarizations)
initial_temp = 300.0;# temperature initial condition
laser_power = 1.0;  # laser power, qdblprime = laser_power*(1+sin(laser_freq*t*two*pi))
spot_size = 4.1e-6; # laser spot size (1/e^2 for Gaussian spot)
laser_freq = 2e7;   # laser modulation frequency
al_thickness = 6e-8;# thickness of transducer layer
dt = 1e-12;         # time step size
nsteps = 100;       # total time steps

use_alsi = true;   # true->ALSI BC on top, false->isothermal BC on top
##############################################################################

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_EXPLICIT)
setSteps(dt, nsteps);

useCUDA(); # This will attempt to use CUDA

# direction and band numbers
ndirs = nphi*ntheta; # total number of directions
(t_bands, l_bands) = get_band_distribution(frequency_bands);
total_bands = t_bands + l_bands; # total band/polarization pairs

# This is needed because we are not partitioning the mesh, only bands.
np = Finch.finch_state.config.num_procs;
rank = Finch.finch_state.config.proc_rank;
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
MPI = Finch.MPI; # This is used in temperature update.

# A simple mesh is internally generated for convenience
# This is a uniform grid of rectangles
mesh(QUADMESH, # quad elements
    elsperdim=[120, 120], # elements in each direction of the uniform grid
    interval=[0, 525e-6, 0, 525e-6],  # interval in each direction
    bids=4, # 4 boundary IDs for this mesh correspond to left, right, bottom, top
    partitions=1) # no partitioning

# area and volume need to be adjusted for axisymmetric
adjust_area_volume();

# Indices, Variables, etc.
direction = index("direction", range = [1,ndirs])
band = index("band", range = [1,nbands])

# These are all set as variables because they are unknown, but only I is solved for in the PDE.
I = variable("I", type=VAR_ARRAY, location=CELL, index = [direction, band]) # Intensity
Io = variable("Io", type=VAR_ARRAY, location=CELL, index = [band]) # Equilibrium intensity for each band
Iplus = variable("Iplus", type=VAR_ARRAY, location=CELL, index = [direction, band]) # A value needed for the axisymmetric case
beta = variable("beta", type=VAR_ARRAY, location=CELL, index = [band]) # Relaxation time scale
temperature = variable("temperature", location=CELL) # temperature of each cell
temperatureLast = variable("temperatureLast", location=CELL) # temperature from last time step
G_last = variable("G_last", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity from last step
G_next = variable("G_next", type=VAR_ARRAY, location=CELL, index = [band]) # integrated intensity for current step

# Coefficients and related numbers
(dir_x, dir_y, dir_z, int_x, int_y, int_z, Omega) = get_directions(nphi, ntheta, true)
(center_freq, delta_freq, polarizations) = get_band_frequencies(frequency_bands, low_band, high_band);
group_v = get_group_speeds(center_freq, polarizations);
(axial_p, axial_pplus) = get_axialp(int_y, nphi, ntheta);

# Pieces needed for ALSI boundary
alsi_info = get_alsi_info(dir_y, dir_z, int_y, int_z, initial_temp);

# These are set as coefficients because they have known values.
# Note that the axisymmetric case uses cylindrical coordinates
Sz = coefficient("Sz", dir_z, type=VAR_ARRAY) # direction z component
Sr = coefficient("Sr", dir_y, type=VAR_ARRAY) # direction radial component
Sp = coefficient("Sp", dir_x, type=VAR_ARRAY) # direction phi component
iSz = coefficient("iSz", int_z, type=VAR_ARRAY) # integrated direction x component
iSr = coefficient("iSr", int_y, type=VAR_ARRAY) # integrated direction y component
omega = coefficient("omega", Omega, type=VAR_ARRAY) # omega
vg = coefficient("vg", group_v, type=VAR_ARRAY) # group speed
axialp = coefficient("axialp", axial_p, type=VAR_ARRAY) # a factor for the axisymmetric case
axialpplus = coefficient("axialpplus", axial_pplus, type=VAR_ARRAY) # a factor for the axisymmetric case

# This extra storage is allocated for parallel temperature update
num_cells = length(temperature.values);
uold = zeros(num_cells); # These strange names are taken from the fortran code
gnb = zeros(num_cells);
unew = zeros(num_cells);
uchange = zeros(num_cells);
gna = zeros(num_cells);
uprime = zeros(num_cells);
converged = fill(false, num_cells);

# bottom symmetric, left and top isothermal, right depends on use_alsi
boundary(I, 1, FLUX, "isothermal_bdry_axi(I, vg, Sz, Sr, iSz, iSr, band, direction, omega, normal, 300)") # left (ID=1)
if use_alsi
    boundary(I, 2, FLUX, "alsi_bdry_axi(I, vg, Sz, Sr, iSz, iSr, band, direction, omega, normal, faceID)") # right (ID=2)
else
    boundary(I, 2, FLUX, "isothermal_bdry_axi(I, vg, Sz, Sr, iSz, iSr, band, direction, omega, normal, 350)") # right (ID=2)
    # boundary(I, 4, FLUX, "isothermal_bdry(I, vg, Sz, Sr, band, direction, normal, (300 + 50*exp(-x*x/(5e-9))))") # right (ID=2)
end
boundary(I, 3, FLUX, "symmetric_bdry_axi(I, vg, Sz, Sr, Sp, iSz, iSr, band, direction, omega, normal)") # bottom (ID=3)
boundary(I, 4, FLUX, "isothermal_bdry_axi(I, vg, Sz, Sr, iSz, iSr, band, direction, omega, normal, 300)") # top (ID=4)


initial(I, [equilibrium_intensity(center_freq[b], delta_freq[b], initial_temp, polarizations[b]) for d=1:ndirs, b=1:nbands])
initial(Io, [equilibrium_intensity(center_freq[b], delta_freq[b], initial_temp, polarizations[b]) for b=1:nbands])
initial(beta, [get_time_scale(center_freq[b], initial_temp, polarizations[b]) for b=1:nbands])
initial(temperature, initial_temp);
initial(temperatureLast, initial_temp);

# assemblyLoops([band, "elements", direction])
assemblyLoops(["elements", band, direction])

# Get integrated intensity from initial equilibrium.
# To get initial values here we have to manually initialize.
evalInitialConditions();
get_integrated_intensity_3d!(G_last.values, I.values, ndirs, nbands, Omega);
update_Iplus!(Iplus.values, I.values, nphi, ntheta, nbands);

# After each time step the temperature, equilibrium I, and time scales are updated
function post_step()
    update_temperature(temperature.values, temperatureLast.values, I.values, beta.values, 
                        center_freq, delta_freq, converged, polarizations, band_parallel=true, threed=true, omega=Omega);
    update_Iplus!(Iplus.values, I.values, nphi, ntheta, nbands);
    if use_alsi
        update_alsi_temp!(alsi_info, I.values, center_freq, delta_freq, polarizations, np);
    end
end
postStepFunction(post_step);

# axisymmetric BTE:
conservationForm(I, "(Io[band] - I[direction,band]) * beta[band] "*
                        "+ vg[band] / (y * omega[direction]) * (axialpplus[direction] * Iplus[direction,band] - axialp[direction] * I[direction,band])"*
                        "+ surface(vg[band] / omega[direction]"*
                        "* upwind([Sz[direction];Sr[direction]] , I[direction,band])"*
                        "* (dot(normal(), [iSz[direction];iSr[direction]])) / (dot(normal(), [Sz[direction]; Sr[direction]])))")


exportCode("bte2daxicode") # uncomment to export generated code to a file
# importCode("bte2daxicodein") # uncomment to import code from a file

solve(I)

##################################################
# End warmup
# If timing after a warmup run, set the steps above to 2
# and uncomment below
##################################################

# if Finch.finch_state.config.proc_rank == 0
#     println("Warm-up result:")
#     show(Finch.finch_state.timer_output)
#     println();
# end
# Finch.reset_timer!(Finch.finch_state.timer_output)

# setSteps(dt, nsteps);

# evalInitialConditions();
# get_integrated_intensity!(G_last.values, I.values, ndirs, nbands);

# solve(I)

# outputValues(temperature, "bte2daxiTemp", format="vtk");

finalizeFinch()

println("max "*string(maximum(temperature.values))*" , min "*string(minimum(temperature.values)))

# display(temperature.values)

# using Plots
# pyplot();
# xy = Finch.finch_state.fv_info.cellCenters
# display(plot(xy[1,:], xy[2,:], temperature.values[:], st=:surface))