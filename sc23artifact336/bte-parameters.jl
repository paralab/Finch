# BTE parameter setup

# Constants #################################################################
# These numbers were adapted from constants_module.f90 from the original code

# Fundamental
const dirac = 1.054571628e-34;
const boltzman = 1.3806503e-23;
const hobol = 7.63822401661014e-12; # dirac/boltzman

# Material
# NOTE: many of these have a 3-letter suffix like "_LAS"
#   1st: L = longitudinal polarization, T = transverse polarization
#   2nd: A = accoustic phonons, O = optical phonons
#   3rd: S = silicon, G = germanium
# For this test we are only considering TAS. Later LAS will be added.
const freq_min_TAS = 0.0;
const freq_max_TAS = 2.97927417405992e13;
const c_TAS = -2.26e-7; # c and vs are coefficients in the wave vector equation:
const vs_TAS= 5230.0;     # wk = wo + vs*|K| + c*|K|^2   where wo=0 for Transverse.

const freq_min_LAS = 0;
const freq_max_LAS = 7.728337675901222e13;
const c_LAS = -2.0e-7;
const vs_LAS= 9010.0;

# Germanium is not used, but we need these to match freq bands with fortran
const freq_max_LAG = 4.41802754820937e13
const freq_max_TAG = 1.50052084481175e13
const use_germanium_bands = false; # set to true to match old fortran version

# Others to be added as needed

# For relaxation time (non-Broido)
const wmax_half = 2.417e13;
const lattice_const = 5.43e-10;
const rho = 2.33e3;
const btu = 5.5e-18;
const bl = 2.0e-24;
const btn =9.3e-13;

# For Broido relaxation time
const debye_temp = 636.0;
const a_n_TAS = 10.9e-20;
const a_u_TAS = 37.8e-47;
const a_n_LAS = 7.10e-20;
const a_u_LAS = 9.51e-47;

# 5-point gaussian quadrature
const g5xi = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664]; # gaussian quadrature points
const g5wi = [0.23692688505618908, 0.47862867049936647, 0.5688888888888889, 0.47862867049936647, 0.23692688505618908]; # gaussian weights

# # 20-point gaussian quadrature
const g20xi = [0.076526521133497, 0.227785851141645, 0.373706088715419, 0.510867001950827, 0.636053680726515,
                0.74633190646015, 0.839116971822218, 0.912234428251325, 0.963971927277913, 0.993128599185094,
                -0.076526521133497, -0.227785851141645, -0.373706088715419, -0.510867001950827, -0.636053680726515,
                -0.74633190646015, -0.839116971822218, -0.912234428251325, -0.963971927277913, -0.993128599185094];
#
const g20wi = [0.152753387130725, 0.149172986472603, 0.142096109318382, 0.131688638449176, 0.118194531961518,
                0.10193011981724, 0.083276741576704, 0.062672048334109, 0.040601429800386, 0.017614007139152,
                0.152753387130725, 0.149172986472603, 0.142096109318382, 0.131688638449176, 0.118194531961518,
                0.10193011981724, 0.083276741576704, 0.062672048334109, 0.040601429800386, 0.017614007139152];
#

# This will be set in get_directions()
# Many things change for the axisymmetric case
axisymmetric = false;

############################################################################

# Direction vectors
# input: number of directions in phi and theta
# output: (x_i, y_i, intx_i, inty_i, omega) direction vector components and integration factors for 2D
# OR
# output: (x_i, y_i, z_i, intx_i, inty_i, intz_i, omega) same for axisymmetric or 3D
function get_directions(nphi, ntheta=1, isAxisymmetric=false)
    global axisymmetric = isAxisymmetric;
    ndir = nphi * ntheta;
    dir_x = zeros(ndir);
    dir_y = zeros(ndir);
    dir_z = zeros(ndir);
    int_x = zeros(ndir);
    int_y = zeros(ndir);
    int_z = zeros(ndir);
    omega = zeros(ndir);
    
    dphi = 2*pi/nphi;
    dtheta = pi/ntheta;
    
    # phi is measured clockwise from the y+ axis
    # theta is measured away from z+ axis
    count = 0
    for j=1:ntheta
        theta = (j-0.5) * dtheta;
        for i=1:nphi
            phi = (i-0.5) * dphi;
            
            count += 1;
            omega[count] = 2 * sin(theta) * sin(dtheta*0.5) * dphi;
            
            # direction vector
            dir_x[count] = sin(theta) * sin(phi);
            dir_y[count] = sin(theta) * cos(phi);
            dir_z[count] = cos(theta);
            
            # integration factor
            int_x[count] = sin(phi) * sin(dphi/2) * (dtheta - cos(2*theta) * sin(dtheta));
            int_y[count] = cos(phi) * sin(dphi/2) * (dtheta - cos(2*theta) * sin(dtheta));
            int_z[count] = dphi / 2 * sin(2*theta) * sin(dtheta);
        end
    end
    
    # Simplified for pure 2D
    if ntheta == 1
        return (dir_x, dir_y, int_x, int_y, omega);
        
    else
        return (dir_x, dir_y, dir_z, int_x, int_y, int_z, omega);
    end
end

# How many bands for each polarization
# T overlaps L, so actual number of bands will be NT + NL
# The bands that overlap will be treated as separate and
# can be done in parallel.
# input: number of frequency bands
# output: number of bands for each polarization
function get_band_distribution(n)
    l_bands = n;
    t_bands = Int(floor(n * freq_max_TAS / freq_max_LAS));
    # If n = 40, t_bands = 15
    
    # Which fortran code version is used?
    if use_germanium_bands
        tg = Int(floor(freq_max_TAG / freq_max_LAS * n));
        ts = Int(floor((freq_max_TAS - freq_max_TAG) / freq_max_LAS * (n-1) + 1));
        lg = Int(floor((freq_max_LAG - freq_max_TAS) / freq_max_LAS * (n-1) + 1));
        ls = Int(floor((freq_max_LAS - freq_max_LAG) / freq_max_LAS * (n-1) + 1));
        ls = tg + ts + lg + ls;
        lg = tg + ts + lg;
        ts = tg + ts;
        
    else
        tg = 0;
        lg = 0;
        ts = Int(floor(freq_max_TAS / freq_max_LAS * n));
        ls = Int(floor((freq_max_LAS - freq_max_TAS) / freq_max_LAS * n + 1));
        ls = ts + ls;
    end
    
    if ls != n
        println("bad band distribution LAS = " * string(ls) * "instead of " * string(n))
        ls = n;
    end
    
    if use_germanium_bands
        return (ts, ls, tg, lg);
    else
        return (ts, ls);
    end
end

# UPDATED for polarizations
# Frequency bands
# input: number of bands
# output: centers of frequency bands, width of bands, polarizations of each band
function get_band_frequencies(n, low=0, high=0)
    # Find all frequency band widths and centers.
    # delta_freq is different for L and T polarizations and materials
    if use_germanium_bands
        #######
        # Due to the fortran version using both polarizations and materials in the band setting
        # we need to replicate it here. It is not ideal, but...
        # order = gt, st, gl, sl
        (st, sl, gt, gl) = get_band_distribution(n);
        
        all_centers = zeros(st+sl);
        all_dw = zeros(st+sl);
        
        for i=1:gt
            all_dw[i] = freq_max_TAG / gt;
        end
        for i=(gt+1):st
            all_dw[i] = (freq_max_TAS - freq_max_TAG) / (st-gt);
        end
        for i=(st+1):gl
            all_dw[i] = (freq_max_LAG - freq_max_TAS) / (gl-st);
        end
        for i=(gl+1):sl
            all_dw[i] = (freq_max_LAS - freq_max_LAG) / (sl-gl);
        end
        
    else # Only silicon
        (st, sl) = get_band_distribution(n);
        
        all_centers = zeros(st+sl);
        all_dw = zeros(st+sl);
        
        for i=1:st
            all_dw[i] = freq_max_TAS / st;
        end
        for i=(st+1):sl
            all_dw[i] = (freq_max_LAS - freq_max_TAS) / (sl-st);
        end
    end
    
    all_centers[1] = all_dw[1] * 0.5;
    for i=2:(st+sl)
        all_centers[i] = all_centers[i-1] + (all_dw[i] + all_dw[i-1]) * 0.5;
    end
    
    if low==0
        low = 1;
    end
    if high==0
        high = st + sl;
    end
    
    # band center freq.
    centers = zeros(high-low+1);
    dw = zeros(high-low+1);
    polarizations = fill("X",high-low+1);
    num_T = max(0, min(st, high) - low + 1);
    num_L = max(0, high - max(st+1,low) + 1);
    first_L = max(1, low-st);
    for i=1:num_T
        centers[i] = all_centers[low+i-1];
        dw[i] = all_dw[low+i-1];
        polarizations[i] = "T";
    end
    for i=1:num_L
        centers[num_T + i] = all_centers[first_L+i-1];
        dw[num_T + i] = all_dw[first_L+i-1];
        polarizations[num_T + i] = "L";
    end
    
    return (centers, dw, polarizations);
end

# Group speed
# input: band center frequencies
# output: group speeds for each band
function get_group_speeds(frequency::Vector{Float64}, polarizations::Vector{String})
    n = length(frequency);
    gv = zeros(n);
    for i=1:n
        if polarizations[i] == "T"
            gv[i] = sqrt(vs_TAS*vs_TAS + 4 * c_TAS * frequency[i]);
        else
            gv[i] = sqrt(vs_LAS*vs_LAS + 4 * c_LAS * frequency[i]);
        end
    end
    return gv;
end

# Scattering time scale from relaxtime.f90 (NOT Broido version)
# input: beta array to update, band center frequencies, temperature
# output: band average inverse time scales (beta = 1/tau)
function get_time_scale!(beta::Matrix{Float64}, freq::Vector{Float64}, temp::Array{Float64}, polarizations::Vector{String})
    n = length(temp); # number of cells
    m = length(freq); # number of bands
    @inbounds begin
    for j=1:n # loop over cells
        for i=1:m # loop over bands
            if polarizations[i] == "T"
                beta[i, j] = btn * (temp[j]^4) * freq[i]
                if freq[i] >= wmax_half
                    beta[i, j] += btu * (freq[i]*freq[i]) / sinh(hobol * freq[i] / temp[j]);
                end
            else # "L"
                beta[i, j] = bl * (temp[j]^3) * (freq[i]*freq[i]);
            end
        end
    end
    end#inbounds
end
# Similar to above, but for a single band and temp
function get_time_scale(freq::Float64, temp::Union{Int, Float64}, polarization::String)
    if polarization == "T"
        beta = btn * (temp*temp*temp*temp) * freq
        if freq >= wmax_half
            beta += btu * (freq*freq) / sinh(hobol * freq / temp);
        end
    else
        beta = bl * (temp*temp*temp) * (freq*freq);
    end

    return beta;
end

# Equilibrium intensity
# Io = integral_dw(dirac*freq/(8*pi^3) * 1/(a) * (b)^2 * dw)
#   where a = exp(hobol*freq/temp) - 1
#         b = (-vs + sqrt(vs^2 + 4*freq*c)) / (2*c)
#
# input: Intensity array to update, band center frequency, band width, temperature of each cell
# output: equilibrium intensity in each cell
function equilibrium_intensity!(intensity::Matrix{Float64}, freq::Vector{Float64}, dw::Vector{Float64}, temp::Array{Float64}, polarizations::Vector{String})
    
    n = length(temp); # should be number of cells
    # dirac/(32*pi^3) = 1.062861036647414e-37
    
    @inbounds begin
    for ci=1:n # loop over cells
        for i=1:length(freq) # loop over bands
            if polarizations[i] == "T"
                vs = vs_TAS;
                c = c_TAS;
                extra_factor = 2;
                const_part = 1.062861036647414e-37 / (c*c) * dw[i]/2; # constants to pull out of integral
            else # "L"
                vs = vs_LAS;
                c = c_LAS;
                extra_factor = 1;
                const_part = 1.062861036647414e-37 / (c*c) * dw[i]/2; # constants to pull out of integral
            end
            tmp = 0.0;
            for gi=1:20
                fi = freq[i] + dw[i]/2 * g20xi[gi]; # frequency at gauss point
                tmp += (fi * (-vs + sqrt(vs*vs + 4*fi*c))^2 / (exp(hobol*fi/temp[ci]) - 1)) * g20wi[gi] * extra_factor;
            end
            intensity[i,ci] = tmp * const_part;
        end
    end
    end#inbounds
end
# Similar to above, but for a single frequency band and temperature
function equilibrium_intensity(freq::Float64, dw::Float64, temp::Union{Int, Float64}, polarization::String)
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
        extra_factor = 2;
    else # L
        vs = vs_LAS;
        c = c_LAS;
        extra_factor = 1;
    end

    # dirac/(32*pi^3) = 1.062861036647414e-37
    const_part = 1.062861036647414e-37 * dw/2 / (c*c); # constants to pull out of integral
    intensity = 0.0;
    @inbounds begin
    for gi=1:20
        fi = freq + dw/2 * g20xi[gi]; # frequency at gauss point
        K2 = (-vs + sqrt(vs*vs + 4*fi*c))^2; # K^2 * (2*c)^2   the (2*c)^2 is put in the const_part
        intensity += (fi * K2 / (exp(hobol*fi/temp) - 1)) * g20wi[gi] * extra_factor;
    end
    end#inbounds
    intensity *= const_part;

    return intensity;
end

# The integrated intensity in each band in each cell (integrated over directions)
function get_integrated_intensity!(int_intensity::Matrix{Float64}, intensity::Matrix{Float64}, ndirs::Int, nbands::Int)
    n = size(intensity,2); # num cells
    omega = 2*pi/ndirs * 2; # angle per direction * 2 
    # Integrate over local bands and cells
    @inbounds begin
    for i=1:n
        for j=1:nbands
            int_intensity[j,i] = 0.0;
            for k=1:ndirs
                int_intensity[j,i] += intensity[(j-1)*ndirs+k, i] * omega;
            end
        end
    end
    end#inbounds
end
function get_integrated_intensity_3d!(int_intensity::Matrix{Float64}, intensity::Matrix{Float64}, ndirs::Int, nbands::Int, omega::Vector{Float64})
    n = size(intensity,2); # num cells
    # Integrate over local bands and cells
    for i=1:n
        for j=1:nbands
            int_intensity[j,i] = 0.0;
            for k=1:ndirs
                int_intensity[j,i] += intensity[(j-1)*ndirs+k, i] * omega[k];
            end
        end
    end
end

# Derivative of equilibrium intensity with respect to temperature
# input: frequency bands, bandwidth, temperature at each cell
# output: dI/dT at each cell
# NOTE: This produces an array of dIdT for every cell and band. Single cell/band version is below
function dIdT(freq::Vector{Float64}, dw::Vector{Float64}, temp::Vector{Float64}, polarizations::Vector{String})
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    n = length(temp); # number of cells
    m = length(freq); # number of bands
    didt = zeros(m, n);
    const_part = 0.5 * dirac * hobol / (8*pi^3);

    @inbounds begin
    for i=1:n # loop over cells
        for j=1:m # loop over bands
            if polarizations[j]=="T"
                vs = vs_TAS;
                c = c_TAS;
                extra_factor = 2;
            else # L
                vs = vs_LAS;
                c = c_LAS;
                extra_factor = 1;
            end
            tmp = 0.0;
            for gi=1:29 # gaussian quadrature
                fi = freq[j] + dw[j]/2 * g20xi[gi]; # frequency at gauss point
                # K2 = ((-vs + sqrt(vs*vs + 4*fi*c)) / (2*c))^2; # K^2
                tmpK = (-vs + sqrt(vs*vs + 4*fi*c)) / (2*c); # K updated to match ipcalc
                tmp2 = exp(hobol*fi/temp[i]);
                # tmp += (tmp2 * fi * K2 / (tmp2 - 1)^2) * wi[gi];
                tmp += (tmp2 * (fi * tmpK)^2 / (tmp2 - 1)^2) * g20wi[gi] * extra_factor; # updated to match ipcalc
            end
            didt[j, i] = const_part * dw[j] * tmp / temp[i]^2;
        end
    end
    end#inbounds

    return didt;
end

# Derivative of equilibrium intensity with respect to temperature
# input: center frequency, bandwidth, temperature at one cell
# output: dI/dT for one cell and band
# NOTE: This is for one cell and band. Full cells/bands version above.
function dIdT_single(freq::Float64, dw::Float64, temp::Float64, polarization::String)
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
        extra_factor = 2;
    else # L
        vs = vs_LAS;
        c = c_LAS;
        extra_factor = 1;
    end

    tmp = 0.0;
    @inbounds begin
    for gi=1:20 # gaussian quadrature
        fi = freq + dw/2 * g20xi[gi]; # frequency at gauss point
        # K2 = ((-vs + sqrt(vs*vs + 4*fi*c)) / (2*c))^2; # K^2
        tmpK = (-vs + sqrt(vs*vs + 4*fi*c)) / (2*c); # K updated to match ipcalc
        tmp2 = exp(hobol*fi/temp);
        # tmp += (tmp2 * fi * K2 / (tmp2 - 1)^2) * wi[gi];
        tmp += (tmp2 * (fi * tmpK)^2 / (tmp2 - 1)^2) * g20wi[gi] * extra_factor; # updated to match ipcalc
    end
    end#inbounds
    didt = tmp * 0.5 * dw * 3.2473482785757725e-48 / (temp * temp); # dirac * hobol / (8*pi^3) = 3.2473482785757725e-48

    return didt;
end

# Change temperature for one time step
# input: temperature array to update, previous step temperature, intensity from this and the previous time step,
#        temperature from previous time step, frequency bands, bandwidth
# output: temperature for next step for each cell
function get_next_temp!(temp_next::Array{Float64}, temp_last::Array{Float64}, I_next::Array{Float64},
                        freq::Vector{Float64}, dw::Vector{Float64},
                        group_v::Vector{Float64}, Io_values::Array{Float64}, G_last_values::Array{Float64}, G_next_values::Array{Float64},
                        polarizations::Vector{String}; threed::Bool=false, omega=nothing)
    debug = false;
    n = length(temp_last); # number of cells
    m = length(freq); # number of bands

    # G_last = get_integrated_intensity(I_last, ndirs, nbands);
    # G_next = get_integrated_intensity(I_next, ndirs, nbands);
    if threed
        get_integrated_intensity_3d!(G_next.values, I_next, ndirs, nbands, omega);
    else
        get_integrated_intensity!(G_next.values, I_next, ndirs, nbands);
    end

    # didt = dIdT(freq, dw, temp_last, polarization=polarization); # Use single version in loop instead

    idt::Float64 = 1.0/dt;

    ave_iters = 0;

    # These tolerances  and iters were hard coded into the fortran desite input parameters
    abstol = 1e-6; # absolute change in U
    contol = 1e-3; # change in U / first change in U
    maxiters = 5;

    @inbounds begin
    for i=1:n # loop over cells
        uold = 0.0;
        uchange = 0.0;
        gnb = 0.0;
        gna = 0.0;
        for j=1:m # loop over bands
            beta = get_time_scale(freq[j], temp_next[i], polarizations[j]);

            uold += Io_values[j,i] / group_v[j];
            uchange += Io_values[j,i] * beta / group_v[j];
            gnb += G_last_values[j,i] * (idt - beta) / group_v[j];
            gna += G_next_values[j,i] * idt / group_v[j];

            G_last_values[j,i] = G_next_values[j,i]; # That was the only place we need G_last, so update it here.
        end
        uold = 4*pi*idt * uold[offset_i];
        uchange = 4*pi * uchange[offset_i];
        
        unew = (uold + gna - gnb - uchange) / idt
        
        delta_T = 0;
        for iter=1:maxiters # iteratively refine delta_T
            uchange = 0.0;
            uprime = 0.0;
            
            for j=1:m # loop over bands
                didt = dIdT_single(freq[j], dw[j], temp_next[i], polarizations[j]);
                
                uchange += equilibrium_intensity(freq[j], dw[j], temp_next[i], polarizations[j]) / group_v[j];
                uprime += didt / group_v[j];
            end
            uchange = 4*pi*uchange;
            uprime = 4*pi*uprime;
            
            delta_u = unew - uchange;
            delta_T = delta_u / uprime;
            
            temp_next[i] = temp_next[i] + delta_T;
            
            # hijack uold for holding uchange_1
            if iter == 1
                uold = delta_u;
            end
            if abs(delta_u) < abstol || abs(delta_u / uold) < contol
                ave_iters += iter;
            else
                if iter==maxiters
                    println("Temperature change didn't converge. Last delta_T = "*string(delta_T));
                end
            end
        end# iterative refinement
    end# cell loop
    end#inbounds
    ave_iters = ave_iters/n;
    if debug
        println("ave iterations: "*string(ave_iters))
    end

    return ave_iters;
end

# Band-based parallel version of above
# Change temperature for one time step
# input: temperature array to update, previous step temperature, intensity from this and the previous time step,
#        temperature from previous time step, frequency bands, bandwidth
# output: temperature for next step for each cell
function get_next_temp_par!(temp_next::Array{Float64}, temp_last::Array{Float64}, I_next::Array{Float64},
                            freq::Vector{Float64}, dw::Vector{Float64},
                            uold::Vector{Float64}, unew::Vector{Float64}, gna::Vector{Float64}, gnb::Vector{Float64}, uprime::Vector{Float64},
                            uchange::Vector{Float64},
                            group_v::Vector{Float64}, Io_values::Array{Float64}, G_last_values::Array{Float64}, G_next_values::Array{Float64},
                            converged::Vector{Bool}, polarizations::Vector{String}; threed::Bool=false, omega=nothing)
    debug = false;
    n = length(temp_last); # number of cells
    m = length(freq); # number of local bands
    nband_partitions::Int = MPI.Comm_size(MPI.COMM_WORLD);

    if threed
        get_integrated_intensity_3d!(G_next_values, I_next, ndirs, nbands, omega);
    else
        get_integrated_intensity!(G_next_values, I_next, ndirs, nbands);
    end
    # didt = dIdT(freq, dw, temp_last, polarization=polarization); # Use single version in loop instead

    idt::Float64 = 1.0/dt;

    ave_iters = 0;

    # These tolerances  and iters were hard coded into the fortran desite input parameters
    abstol = 1e-6; # absolute change in U
    contol = 1e-3; # change in U / first change in U
    maxiters = 20;

    # First find uold and gnb from the previous step
    rank = MPI.Comm_rank(MPI.COMM_WORLD);
    @inbounds for i=1:n # loop over cells
        uold[i] = 0.0;
        uchange[i] = 0.0;
        gnb[i] = 0.0;
        gna[i] = 0.0;
        for j=1:m # loop over local bands
            beta = get_time_scale(freq[j], temp_next[i], polarizations[j]);

            uold[i] += Io_values[j,i] / group_v[j];
            uchange[i] += Io_values[j,i] * beta / group_v[j];
            gnb[i] += G_last_values[j,i] * (idt - beta) / group_v[j];
            gna[i] += G_next_values[j,i] * idt / group_v[j];
            
            G_last_values[j,i] = G_next_values[j,i]; # That was the only place we need G_last, so update it here.
            
        end
        uold[i] = 4*pi*idt * uold[i];
        uchange[i] = 4*pi * uchange[i];

        unew[i] = (uold[i] + gna[i] - gnb[i] - uchange[i]) / idt
    end
    
    # Now unew holds the partially reduced values for this proc.
    # Use MPI.Allreduce! to do the full reduction and send to all procs.
    MPI.Allreduce!(unew, MPI.SUM, MPI.COMM_WORLD);
    
    # Then iteratively refine delta_T
    converged .= false;
    num_not_converged = n;
    for iter=1:maxiters
        # The cell loop is inside the iteration loop so that gathers can be done for all cells rather than each cell
        @inbounds for i=1:n # loop over cells
            if converged[i]
                continue;
            end
            uchange[i] = 0.0;
            uprime[i] = 0.0;

            for j=1:m # loop over local bands
                beta = get_time_scale(freq[j], temp_next[i], polarizations[j]);
                didt = dIdT_single(freq[j], dw[j], temp_next[i], polarizations[j]);

                uchange[i] += equilibrium_intensity(freq[j], dw[j], temp_next[i], polarizations[j]) / group_v[j];
                uprime[i] += didt / group_v[j];
            end
            uchange[i] = 4*pi*uchange[i];
            uprime[i] = 4*pi*uprime[i];
        end
        
        MPI.Allreduce!(uchange, MPI.SUM, MPI.COMM_WORLD);
        MPI.Allreduce!(uprime, MPI.SUM, MPI.COMM_WORLD);

        # update temp
        @inbounds for i=1:n # loop over cells
            if converged[i]
                continue;
            end

            delta_u = unew[i] - uchange[i];
            delta_T = delta_u / uprime[i];

            temp_next[i] = temp_next[i] + delta_T;
            
            # hijack uold for holding uchange_1
            if iter == 1
                uold[i] = delta_u;
            end
            if abs(delta_u) < abstol || abs(delta_u / uold[i]) < contol
                ave_iters += iter;
                converged[i] = true;
                num_not_converged = num_not_converged - 1;
            else
                if iter==maxiters
                    println("Temperature change didn't converge. Last delta_T = "*string(delta_T));
                end
            end
        end# cell loop

        if num_not_converged == 0
            break;
        end
    end# iterative refinement
    ave_iters = ave_iters/n;
    if debug
        println("ave iterations: "*string(ave_iters))
    end

    return ave_iters;
end

# Update temperature, equilibrium I, and time scale
#
function update_temperature(temp::Matrix{Float64}, temp_last::Matrix{Float64}, I_next::Matrix{Float64}, beta::Matrix{Float64},
                                freq::Vector{Float64}, dw::Vector{Float64}, converged::Vector{Bool}, polarizations::Vector{String}; 
                                band_parallel=false, threed=false, omega=nothing)
    ncells = length(temp);
    @inbounds for i=1:ncells
        temp_last[i] = temp[i];
    end

    if band_parallel
        iterations = get_next_temp_par!(temp, temp_last, I_next, freq, dw, uold, unew, gna, gnb, uprime, uchange,
                                        group_v, Io.values, G_last.values, G_next.values, converged, polarizations, threed=threed, omega=omega);
    else
        iterations = get_next_temp!(temp, temp_last, I_next, freq, dw,
                                    group_v, Io.values, G_last.values, G_next.values, polarizations, threed=threed, omega=omega);
    end

    equilibrium_intensity!(Io.values, freq, dw, temp, polarizations);
    get_time_scale!(beta, freq, temp, polarizations);

    # check for problems
    # for i=1:length(I_next)
    #     if I_next[i] === NaN
    #         println("Error: NaN found in I");
    #         exit(0);
    #     end
    # end

    return iterations;
end

# This is a quantity needed for the axisymmetric case
function get_axialp(int_y, nphi, ntheta)
    ndirs = nphi * ntheta;
    axialp = zeros(ndirs);
    
    count = 0
    for j=1:ntheta
        for i=1:nphi
            count += 1;
            if i == 1
                axialp[count] = 0.0;
            else
                axialp[count] = axialp[count-1] + int_y[count-1];
            end
        end
    end
    
    return axialp;
end

# This is needed for the axisymmetric case
function update_Iplus!(Iplus, I, nphi, ntheta, nbands)
    ncells = size(I,2);
    
    for celli = 1:ncells
        count = 0
        for bandi = 1:nbands
            for j=1:ntheta
                for i=1:nphi
                    count += 1;
                    if i == nphi
                        Iplus[count, celli] = 0.0;
                    else
                        Iplus[count, celli] = I[count + 1, celli];
                    end
                end
            end
        end
    end
end