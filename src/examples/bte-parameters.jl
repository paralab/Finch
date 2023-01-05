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

const freq_min_LAS = freq_max_TAS; # In the original code this was 0??
const freq_max_LAS = 7.728337675901222e13;
const c_LAS = -2.0e-7;
const vs_LAS= 9010.0;

const polarization_ratio = 1;

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
# g20xi = zeros(20);
# g20wi = zeros(20);
# g20xi[1] = 0.076526521133497; g20wi[1] = 0.152753387130725;
# g20xi[2] = 0.227785851141645; g20wi[2] = 0.149172986472603;
# g20xi[3] = 0.373706088715419; g20wi[3] = 0.142096109318382;
# g20xi[4] = 0.510867001950827; g20wi[4] = 0.131688638449176;
# g20xi[5] = 0.636053680726515; g20wi[5] = 0.118194531961518;
# g20xi[6] = 0.746331906460150; g20wi[6] = 0.101930119817240;
# g20xi[7] = 0.839116971822218; g20wi[7] = 0.083276741576704;
# g20xi[8] = 0.912234428251325; g20wi[8] = 0.062672048334109;
# g20xi[9] = 0.963971927277913; g20wi[9] = 0.040601429800386;
# g20xi[10] = 0.993128599185094; g20wi[10] = 0.017614007139152;
# for i = 1:10
#     g20xi[10+i] = -g20xi[i];
#     g20wi[10+i] = g20wi[i];
# end
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

############################################################################

# Direction vectors
# input: number of directions
# output: (x_i, y_i) direction vector components
function get_directions_2d(n)
    dir_x = zeros(n);
    dir_y = zeros(n);
    reflect = zeros(Int,n); # The index of the reflecting direction
    dphi = 2*pi/n;
    phione = dphi/2;
    for di=1:n
        # dir_x[di] = cos(2*pi*(di-1)/n);
        # dir_y[di] = sin(2*pi*(di-1)/n);

        # modified to match original
        thisphi = phione + (di-1)*dphi;
        dir_x[di] = sin(thisphi);
        dir_y[di] = cos(thisphi);
    end
    if mod(n,2) == 0
        halfn = Int(n/2);
        for di=1:halfn
            reflect[di] = n-di+1;
            reflect[n-di+1] = di;
        end
    else
        println("Please use an even number of directions for reflection purposes. Reflections will be incorrect.")
        halfn = Int(floor(n/2));
        for di=1:halfn
            reflect[di] = n-di+1;
            reflect[n-di+1] = di;
        end
        reflect[halfn+1] = 1;
    end

    return (dir_x, dir_y, reflect);
end

# Direction vectors
# input: number of directions
# output: (x_i, y_i, z_i) direction vector components
#         omega
function get_directions_3d(nphi, ntheta)
    n = nphi * ntheta;
    dir_x = zeros(n);
    dir_y = zeros(n);
    dir_z = zeros(n);
    omega = zeros(n);

    dphi = 2*pi/nphi;
    dtheta = pi/ntheta;
    ind = 1;
    for i=1:ntheta
        theta = (i-0.5)*dtheta;
        for j=1:nphi
            anphi = (j-0.5)*dphi;
            dir_x[ind] = sin(theta)*sin(anphi);
            dir_y[ind] = sin(theta)*cos(anphi);
            dir_z[ind] = cos(theta);
            omega[ind] = 2*sin(theta)*sin(0.5*dtheta)*dphi;
            # inrmu[ind] = sin(anphi) * sin(dphi*0.5) * (dtheta - cos(2*theta)*sin(dtheta));
            # inrxi[ind] = cos(anphi) * sin(dphi*0.5) * (dtheta - cos(2*theta)*sin(dtheta));
            # inret[ind] = 0.5 * dphi * sin(2*theta) * sin(dtheta);
            ind += 1;
        end
    end
    return (dir_x, dir_y, dir_z, omega);
end

# For 3d the reflected direction is not so simple to store.
# This computes it.
function get_3d_reflection(dir_x, dir_y, dir_z, S, N)
    sout = 0;
    sdotn = S[1]*N[1] + S[2]*N[2] + S[3]*N[3];

    # ideal reflected direction
    srx = S[1] - 2 * sdotn * N[1];
    sry = S[2] - 2 * sdotn * N[2];
    srz = S[3] - 2 * sdotn * N[3];

    # Find the direction that is closest to the reflected direction
    min_off = 1e6;
    for i = 1:length(dir_x)
        sdots0 = dir_x[i]*srx + dir_y[i]*sry + dir_z[i]*srz;
        if abs(sdots0 - 1) < min_off
            sout = i;
            min_off = abs(sdots0 - 1);
        end
    end

    return sout;
end

# Finds all of the reflection indices for each face of a cube
function get_3d_cube_reflections(dir_x, dir_y, dir_z)
    ndirs = length(dir_x);
    reflect = zeros(Int, ndirs, 6);

    for i=1:ndirs
        S = [dir_x[i], dir_y[i], dir_z[i]];
        # x=0 face
        N = [-1,0,0];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
        # x=1 face
        N = [1,0,0];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
        # y=0 face
        N = [0,-1,0];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
        # y=1 face
        N = [0,1,0];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
        # z=0 face
        N = [0,0,-1];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
        # z=1 face
        N = [0,0,1];
        reflect[i, 1] = get_3d_reflection(dir_x, dir_y, dir_z, S, N);
    end

    return reflect;
end

# Frequency bands
# input: number of bands
# output: centers of frequency bands, width of bands
function get_band_frequencies(n, low=0, high=0; polarization="T")
    if polarization=="T"
        wmin = freq_min_TAS;
        wmax = freq_max_TAS;
    else # L
        wmin = freq_min_LAS;
        wmax = freq_max_LAS;
    end

    if low==0
        low = 1;
    end
    if high==0
        high = n;
    end

    dw = (wmax-wmin)/n;
    centers = zeros(high-low+1);
    for i=1:length(centers)
        centers[i] = wmin + (i - 1)*dw + (low-1)*dw + dw/2;
    end
    return (centers, dw);
end

# Group speed
# input: band center frequencies
# output: group speeds for each band
function get_group_speeds(freq; polarization="T")
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
    else # L
        vs = vs_LAS;
        c = c_LAS;
    end

    n = length(freq);
    gv = zeros(n);
    for i=1:n
        gv[i] = sqrt(vs*vs + 4 * c * freq[i]);
    end
    return gv;
end

# Scattering time scale from relaxtime.f90 (NOT Broido version)
# input: beta array to update, band center frequencies, temperature
# output: band average inverse time scales (beta = 1/tau)
function get_time_scale!(beta::Matrix{Float64}, freq::Vector{Float64}, temp::Array{Float64}; polarization="T")
    n = length(temp); # number of cells
    m = length(freq); # number of bands
    @inbounds begin
    for j=1:n # loop over cells
        for i=1:m # loop over bands
            if polarization == "T"
                beta[i, j] = btn * (temp[j]^4) * freq[i]
                if freq[i] >= wmax_half
                    beta[i, j] += btu * (freq[i]*freq[i]) / sinh(hobol * freq[i] / temp[j]);
                end
            else
                beta[i, j] = bl * (temp[j]^3) * (freq[i]*freq[i]);
            end
        end
    end
    end#inbounds
end
# Similar to above, but for a single band and temp
function get_time_scale(freq::Float64, temp::Union{Int, Float64}; polarization="T")
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
function equilibrium_intensity!(intensity::Matrix{Float64}, freq::Vector{Float64}, dw::Float64, temp::Array{Float64}; polarization="T")
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

    n = length(temp); # should be number of cells
    # dirac/(32*pi^3) = 1.062861036647414e-37
    const_part = 1.062861036647414e-37 / (c*c) * dw/2; # constants to pull out of integral
    @inbounds begin
    for ci=1:n # loop over cells
        for i=1:length(freq) # loop over bands
            tmp = 0.0;
            for gi=1:20
                fi = freq[i] + dw/2 * g20xi[gi]; # frequency at gauss point
                tmp += (fi * (-vs + sqrt(vs*vs + 4*fi*c))^2 / (exp(hobol*fi/temp[ci]) - 1)) * g20wi[gi] * extra_factor;
            end
            intensity[i,ci] = tmp * const_part;
        end
    end
    end#inbounds
end
# Similar to above, but for a single frequency band and temperature
function equilibrium_intensity(freq::Float64, dw::Float64, temp::Union{Int, Float64}; polarization="T")
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
function get_integrated_intensity_3d!(int_intensity, intensity, ndirs, nbands, omega)
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
function dIdT(freq, dw, temp; polarization="T")
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
    else # L
        vs = vs_LAS;
        c = c_LAS;
    end

    n = length(temp); # number of cells
    m = length(freq); # number of bands
    didt = zeros(m, n);
    const_part = dirac * hobol / (8*pi^3);

    @inbounds begin
    for i=1:n # loop over cells
        for j=1:m # loop over bands
            tmp = 0.0;
            for gi=1:29 # gaussian quadrature
                fi = freq[j] + dw/2 * g20xi[gi]; # frequency at gauss point
                # K2 = ((-vs + sqrt(vs*vs + 4*fi*c)) / (2*c))^2; # K^2
                tmpK = (-vs + sqrt(vs*vs + 4*fi*c)) / (2*c); # K updated to match ipcalc
                tmp2 = exp(hobol*fi/temp[i]);
                # tmp += (tmp2 * fi * K2 / (tmp2 - 1)^2) * wi[gi];
                tmp += (tmp2 * (fi * tmpK)^2 / (tmp2 - 1)^2) * g20wi[gi]; # updated to match ipcalc
            end
            didt[j, i] = const_part * dw * tmp / temp[i]^2;
        end
    end
    end#inbounds

    return didt;
end

# Derivative of equilibrium intensity with respect to temperature
# input: center frequency, bandwidth, temperature at one cell
# output: dI/dT for one cell and band
# NOTE: This is for one cell and band. Full cells/bands version above.
function dIdT_single(freq::Float64, dw::Float64, temp::Float64; polarization="T")
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
    else # L
        vs = vs_LAS;
        c = c_LAS;
    end

    tmp = 0.0;
    @inbounds begin
    for gi=1:20 # gaussian quadrature
        fi = freq + dw/2 * g20xi[gi]; # frequency at gauss point
        # K2 = ((-vs + sqrt(vs*vs + 4*fi*c)) / (2*c))^2; # K^2
        tmpK = (-vs + sqrt(vs*vs + 4*fi*c)) / (2*c); # K updated to match ipcalc
        tmp2 = exp(hobol*fi/temp);
        # tmp += (tmp2 * fi * K2 / (tmp2 - 1)^2) * wi[gi];
        tmp += (tmp2 * (fi * tmpK)^2 / (tmp2 - 1)^2) * g20wi[gi]; # updated to match ipcalc
    end
    end#inbounds
    didt = tmp * dw * 3.2473482785757725e-48 / (temp * temp); # dirac * hobol / (8*pi^3) = 3.2473482785757725e-48

    return didt;
end

# Change temperature for one time step
# input: temperature array to update, previous step temperature, intensity from this and the previous time step,
#        temperature from previous time step, frequency bands, bandwidth
# output: temperature for next step for each cell
function get_next_temp!(temp_next::Array{Float64}, temp_last::Array{Float64}, I_next::Array{Float64},
                        freq::Vector{Float64}, dw::Float64,
                        group_v::Vector{Float64}, Io_values::Array{Float64}, G_last_values::Array{Float64}, G_next_values::Array{Float64};
                        polarization::String="T", threed::Bool=false, omega=nothing)
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
            beta = get_time_scale(freq[j], temp_next[i], polarization=polarization);

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
                didt = dIdT_single(freq[j], dw, temp_last[i], polarization=polarization);

                uchange += equilibrium_intensity(freq[j], dw, temp_next[i], polarization=polarization) / group_v[j];
                uprime += 4*pi * didt / group_v[j];
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
                            freq::Vector{Float64}, dw::Float64,
                            uold::Vector{Float64}, unew::Vector{Float64}, gna::Vector{Float64}, gnb::Vector{Float64}, uprime::Vector{Float64},
                            uchange::Vector{Float64},
                            group_v::Vector{Float64}, Io_values::Array{Float64}, G_last_values::Array{Float64}, G_next_values::Array{Float64},
                            converged::Vector{Bool};
                            polarization::String="T", threed::Bool=false, omega=nothing)
    debug = false;
    n = length(temp_last); # number of cells
    m = length(freq); # number of local bands
    nband_partitions::Int = MPI.Comm_size(cell_comm);

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
    maxiters = 5;

    # First find uold and gnb from the previous step
    cell_comm_rank = MPI.Comm_rank(cell_comm);
    offset::Int = cell_comm_rank * n;
    @inbounds for i=1:n # loop over cells
        offset_i = offset + i;
        uold[offset_i] = 0.0;
        uchange[offset_i] = 0.0;
        gnb[offset_i] = 0.0;
        gna[offset_i] = 0.0;
        for j=1:m # loop over local bands
            beta = get_time_scale(freq[j], temp_next[i], polarization=polarization);

            uold[offset_i] += Io_values[j,i] / group_v[j];
            uchange[offset_i] += Io_values[j,i] * beta / group_v[j];
            gnb[offset_i] += G_last_values[j,i] * (idt - beta) / group_v[j];
            gna[offset_i] += G_next_values[j,i] * idt / group_v[j];

            G_last_values[j,i] = G_next_values[j,i]; # That was the only place we need G_last, so update it here.
        end
        uold[offset_i] = 4*pi*idt * uold[offset_i];
        uchange[offset_i] = 4*pi * uchange[offset_i];

        unew[offset_i] = (uold[offset_i] + gna[offset_i] - gnb[offset_i] - uchange[offset_i]) / idt
    end

    # Allgather across this mesh communicator to combine bands
    MPI.Allgather!(unew_buffer, cell_comm);

    # reduce
    @inbounds for i=1:n # loop over cells
        for j=2:nband_partitions # loop over band partitions
            unew[i] += unew[i + (j-1)*n];
        end
    end

    # Then iteratively refine delta_T
    converged .= false;
    num_not_converged = n;
    for iter=1:maxiters
        # The cell loop is inside the iteration loop so that gathers can be done for all cells rather than each cell
        @inbounds for i=1:n # loop over cells
            if converged[i]
                continue;
            end
            offset_i = offset + i;
            uchange[offset_i] = 0.0;
            uprime[offset_i] = 0.0;

            for j=1:m # loop over local bands
                beta = get_time_scale(freq[j], temp_next[i], polarization=polarization);
                didt = dIdT_single(freq[j], dw, temp_last[i], polarization=polarization);

                uchange[offset_i] += equilibrium_intensity(freq[j], dw, temp_next[i], polarization=polarization) / group_v[j];
                uprime[offset_i] += 4*pi * didt / group_v[j];
            end
            uchange[offset_i] = 4*pi*uchange[offset_i];
            uprime[offset_i] = 4*pi*uprime[offset_i];
        end

        # Allgather across this mesh communicator to combine bands
        MPI.Allgather!(uchange_buffer, cell_comm);
        MPI.Allgather!(uprime_buffer, cell_comm);

        # update temp
        @inbounds for i=1:n # loop over cells
            if converged[i]
                continue;
            end
            # reduce
            for j=2:nband_partitions # loop over band partitions
                uchange[i] += uchange[i + (j-1)*n];
                uprime[i] += uprime[i + (j-1)*n];
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
                                freq::Vector{Float64}, dw::Float64, converged::Vector{Bool}; band_parallel=false, threed=false, omega=nothing)
    ncells = length(temp);
    @inbounds for i=1:ncells
        temp_last[i] = temp[i];
    end

    if band_parallel
        iterations = get_next_temp_par!(temp, temp_last, I_next, freq, dw, uold, unew, gna, gnb, uprime, uchange,
                                        group_v, Io.values, G_last.values, G_next.values, converged, threed=threed, omega=omega);
    else
        iterations = get_next_temp!(temp, temp_last, I_next, freq, dw,
                                    group_v, Io.values, G_last.values, G_next.values, threed=threed, omega=omega);
    end

    equilibrium_intensity!(Io.values, freq, dw, temp);
    get_time_scale!(beta, freq, temp);

    # check for problems
    # for i=1:length(I_next)
    #     if I_next[i] === NaN
    #         println("Error: NaN found in I");
    #         exit(0);
    #     end
    # end

    return iterations;
end

#############################################################################
## alternative models.

# What was this??????

# # Change temperature for one time step
# # input: temperature array to update, previous step temperature, intensity from this and the previous time step,
# #        temperature from previous time step, frequency bands, bandwidth
# # output: temperature for next step for each cell
# function get_next_temp!(temp_next, temp_last, I_next, freq, dw; polarization="T", threed=false, omega=nothing)
#     debug = false;
#     n = length(temp_last); # number of cells
#     m = length(freq); # number of bands

#     # G_last = get_integrated_intensity(I_last, ndirs, nbands);
#     # G_next = get_integrated_intensity(I_next, ndirs, nbands);
#     if threed
#         get_integrated_intensity_3d!(G_next.values, I_next, ndirs, nbands, omega);
#     else
#         get_integrated_intensity!(G_next.values, I_next, ndirs, nbands);
#     end

#     # didt = dIdT(freq, dw, temp_last, polarization=polarization); # Use single version in loop instead

#     groupV::Vector{Float64} = group_v;

#     idt = 1/dt;

#     max_iters = 0;
#     ave_iters = 0;

#     tol = 1e-5;
#     maxiters = 200;
#     @inbounds begin
#     for i=1:n # loop over cells
#         # old values are not updated
#         uold = 0.0; # These strange names are taken from the fortran code
#         gnb = 0.0;  #
#         for j=1:m # loop over bands
#             uold += Io.values[j,i] * idt / groupV[j];
#             gnb += G_last.values[j,i] * idt / groupV[j];
#             G_last.values[j,i] = G_next.values[j,i]; # That was the only place we need G_last, so update it here.
#         end
#         uold = 4*pi*uold;

#         delta_T = 0;
#         for iter=1:maxiters # iteratively refine delta_T
#             unew = 0.0; # These strange names are taken from the fortran code
#             gna = 0.0;  #
#             uprime = 0.0;
#             for j=1:m # loop over bands
#                 beta = 1 / get_time_scale(freq[j], temp_next[i], polarization=polarization);
#                 didt = dIdT_single(freq[j], dw, temp_last[i], polarization=polarization);

#                 unew += equilibrium_intensity(freq[j], dw, temp_next[i], polarization=polarization) * (beta + idt) / groupV[j];
#                 gna += G_next.values[j,i] * (beta + idt) / groupV[j];
#                 uprime += 4*pi * didt * (beta + idt) / groupV[j];
#             end
#             unew = 4*pi*unew;

#             delta_T = (uold + gna - gnb - unew) / uprime;

#             # if debug println("cell "*string(i)*" ("*string(uold)*", "*string(unew)*", "*string(gna)*", "*string(gnb)*", "*string(uprime)*") : "*string(delta_T)) end

#             temp_next[i] = temp_next[i] + delta_T;

#             rel_delta = delta_T / temp_next[i];
#             if abs(rel_delta) < tol
#                 max_iters = max(max_iters, iter);
#                 ave_iters += iter;
#                 break;
#             else
#                 if iter==maxiters
#                     println("Temperature change didn't converge. Last delta_T = "*string(delta_T));
#                 end
#             end
#         end# iterative refinement
#     end# cell loop
#     end#inbounds
#     ave_iters = ave_iters/n;
#     if debug
#         println("ave iterations: "*string(ave_iters)*"  max: "*string(max_iters))
#     end

#     return ave_iters;
# end

# # Band-based parallel version of above
# # Change temperature for one time step
# # input: temperature array to update, previous step temperature, intensity from this and the previous time step,
# #        temperature from previous time step, frequency bands, bandwidth
# # output: temperature for next step for each cell
# function get_next_temp_par!(temp_next::Array{Float64}, temp_last::Array{Float64}, I_next::Array{Float64},
#                             freq::Vector{Float64}, dw::Float64,
#                             uold::Vector{Float64}, unew::Vector{Float64}, gna::Vector{Float64}, gnb::Vector{Float64}, uprime::Vector{Float64},
#                             group_v::Vector{Float64}, Io_values::Array{Float64}, G_last_values::Array{Float64}, G_next_values::Array{Float64};
#                             polarization::String="T", threed::Bool=false, omega=nothing)
#     debug = false;
#     n = length(temp_last); # number of cells
#     m = length(freq); # number of local bands
#     nband_partitions::Int = MPI.Comm_size(cell_comm);

#     if threed
#         get_integrated_intensity_3d!(G_next_values, I_next, ndirs, nbands, omega);
#     else
#         get_integrated_intensity!(G_next_values, I_next, ndirs, nbands);
#     end
#     # didt = dIdT(freq, dw, temp_last, polarization=polarization); # Use single version in loop instead

#     idt::Float64 = 1.0/dt;

#     max_iters = 0;
#     ave_iters = 0;

#     tol = 1e-5; # relative tolerance
#     maxiters = 200;

#     # First find uold and gnb from the previous step
#     cell_comm_rank = MPI.Comm_rank(cell_comm);
#     offset::Int = cell_comm_rank * n;
#     @inbounds for i=1:n # loop over cells
#         offset_i = offset + i;
#         uold[offset_i] = 0.0;
#         gnb[offset_i] = 0.0;
#         for j=1:m # loop over local bands
#             uold[offset_i] += Io_values[j,i] * idt / group_v[j];
#             gnb[offset_i] += G_last_values[j,i] * idt / group_v[j];
#             G_last_values[j,i] = G_next_values[j,i]; # That was the only place we need G_last, so update it here.
#         end
#         uold[offset_i] = 4*pi*uold[offset_i];
#     end

#     # Allgather across this mesh communicator to combine bands
#     MPI.Allgather!(uold_buffer, cell_comm);
#     MPI.Allgather!(gnb_buffer, cell_comm);

#     # reduce
#     @inbounds for i=1:n # loop over cells
#         for j=2:nband_partitions # loop over band partitions
#             uold[i] += uold[i + (j-1)*n];
#             gnb[i] += gnb[i + (j-1)*n];
#         end
#     end

#     # Then iteratively refine delta_T
#     converged = fill(false, n);
#     num_not_converged = n;
#     for iter=1:maxiters
#         @inbounds for i=1:n # loop over cells
#             if converged[i]
#                 continue;
#             end
#             offset_i = offset + i;
#             unew[offset_i] = 0.0;
#             gna[offset_i] = 0.0;
#             uprime[offset_i] = 0.0;

#             for j=1:m # loop over local bands
#                 beta = 1 / get_time_scale(freq[j], temp_next[i], polarization=polarization);
#                 didt = dIdT_single(freq[j], dw, temp_last[i], polarization=polarization);

#                 unew[offset_i] += equilibrium_intensity(freq[j], dw, temp_next[i], polarization=polarization) * (beta + idt) / group_v[j];
#                 gna[offset_i] += G_next_values[j,i] * (beta + idt) / group_v[j];
#                 uprime[offset_i] += 4*pi * didt * (beta + idt) / group_v[j];
#             end
#             unew[offset_i] = 4*pi*unew[offset_i];
#         end

#         # Allgather across this mesh communicator to combine bands
#         MPI.Allgather!(unew_buffer, cell_comm);
#         MPI.Allgather!(gna_buffer, cell_comm);
#         MPI.Allgather!(uprime_buffer, cell_comm);

#         # reduce
#         @inbounds for i=1:n # loop over cells
#             if converged[i]
#                 continue;
#             end
#             for j=2:nband_partitions # loop over band partitions
#                 unew[i] += unew[i + (j-1)*n];
#                 gna[i] += gna[i + (j-1)*n];
#                 uprime[i] += uprime[i + (j-1)*n];
#             end

#             delta_T = (uold[i] + gna[i] - gnb[i] - unew[i]) / uprime[i];

#             temp_next[i] = temp_next[i] + delta_T;

#             rel_delta = delta_T / temp_next[i];
#             if abs(rel_delta) < tol
#                 max_iters = max(max_iters, iter);
#                 ave_iters += iter;
#                 converged[i] = true;
#                 num_not_converged = num_not_converged - 1;
#             else
#                 if iter==maxiters
#                     println("Temperature change didn't converge. Last delta_T = "*string(delta_T));
#                 end
#             end
#         end# cell loop

#         if num_not_converged == 0
#             break;
#         end
#     end# iterative refinement
#     ave_iters = ave_iters/n;
#     if debug
#         println("ave iterations: "*string(ave_iters)*"  max: "*string(max_iters))
#     end

#     return ave_iters;
# end

# # Scattering time scale using Broido data from relaxtime_broido.f90
# # input: tau array to update, band center frequencies, temperature
# # output: band average time scales
# function get_time_scale!(tau::Array, freq::Array, temp::Array; polarization="T")
#     if polarization=="T"
#         a_n = a_n_TAS;
#         a_u = a_u_TAS;
#     else # L
#         a_n = a_n_LAS;
#         a_u = a_u_LAS;
#     end

#     n = length(temp); # number of cells
#     m = length(freq); # number of bands
#     for j=1:n # loop over cells
#         for i=1:m # loop over bands
#             tmp = temp[j]*(1 - exp(-3*temp[j]/debye_temp));
#             tau[i, j] = 1 / (a_n * freq[i]^2 * tmp + a_u * freq[i]^4 * tmp);
#         end
#     end
# end
# # Similar to above, but for a single band and temp
# function get_time_scale(freq::Number, temp::Number; polarization="T")
#     if polarization=="T"
#         a_n = a_n_TAS;
#         a_u = a_u_TAS;
#     else # L
#         a_n = a_n_LAS;
#         a_u = a_u_LAS;
#     end

#     tmp = temp*(1 - exp(-3*temp/debye_temp));
#     tau = 1 / (a_n * freq^2 * tmp + a_u * freq^4 * tmp);

#     return tau;
# end