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
const freq_min_TAS = 0;
const freq_max_TAS = 2.97927417405992e13;
const c_TAS = -2.26e-7; # c and vs are coefficients in the wave vector equation:
const vs_TAS= 5230;     # wk = wo + vs*|K| + c*|K|^2   where wo=0 for Transverse.

const freq_min_LAS = freq_max_TAS; # In the original code this was 0??
const freq_max_LAS = 7.728337675901222e13;
const c_LAS = -2.0e-7;
const vs_LAS= 9010;

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
g5xi = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664]; # gaussian quadrature points
g5wi = [0.23692688505618908, 0.47862867049936647, 0.5688888888888889, 0.47862867049936647, 0.23692688505618908]; # gaussian weights

# 20-point gaussian quadrature
g20xi = zeros(20);
g20wi = zeros(20);
g20xi[1] = 0.076526521133497; g20wi[1] = 0.152753387130725;
g20xi[2] = 0.227785851141645; g20wi[2] = 0.149172986472603;
g20xi[3] = 0.373706088715419; g20wi[3] = 0.142096109318382;
g20xi[4] = 0.510867001950827; g20wi[4] = 0.131688638449176;
g20xi[5] = 0.636053680726515; g20wi[5] = 0.118194531961518;
g20xi[6] = 0.746331906460150; g20wi[6] = 0.101930119817240;
g20xi[7] = 0.839116971822218; g20wi[7] = 0.083276741576704;
g20xi[8] = 0.912234428251325; g20wi[8] = 0.062672048334109;
g20xi[9] = 0.963971927277913; g20wi[9] = 0.040601429800386;
g20xi[10] = 0.993128599185094; g20wi[10] = 0.017614007139152;
for i = 1:10
    g20xi[10+i] = -g20xi[i];
    g20wi[10+i] = g20wi[i];
end


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

# Frequency bands
# input: number of bands
# output: centers of frequency bands, width of bands
function get_band_frequencies(n; polarization="T")
    if polarization=="T"
        wmin = freq_min_TAS;
        wmax = freq_max_TAS;
    else # L
        wmin = freq_min_LAS;
        wmax = freq_max_LAS;
    end
    
    dw = (wmax-wmin)/n;
    centers = zeros(n);
    for i=1:n
        centers[i] = wmin + (n-1)*dw + dw/2;
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
function get_time_scale!(beta::Array, freq::Array, temp::Array; polarization="T")
    n = length(temp); # number of cells
    m = length(freq); # number of bands
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
end
# Similar to above, but for a single band and temp
function get_time_scale(freq::Number, temp::Number; polarization="T")
    if polarization == "T"
        beta = btn * (temp^4) * freq
        if freq >= wmax_half                  
            beta += btu * (freq*freq) / sinh(hobol * freq / temp);             
        end
    else
        beta = bl * (temp^3) * (freq*freq);  
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
function equilibrium_intensity!(intensity::Array, freq::Array, dw, temp::Array; polarization="T")
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
end
# Similar to above, but for a single frequency band and temperature
function equilibrium_intensity(freq::Number, dw, temp::Number; polarization="T")
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
    for gi=1:20
        fi = freq + dw/2 * g20xi[gi]; # frequency at gauss point
        K2 = (-vs + sqrt(vs*vs + 4*fi*c))^2; # K^2 * (2*c)^2   the (2*c)^2 is put in the const_part
        intensity += (fi * K2 / (exp(hobol*fi/temp) - 1)) * g20wi[gi] * extra_factor;
    end
    intensity *= const_part;
    
    return intensity;
end

# The integrated intensity in each band in each cell (integrated over directions)
function get_integrated_intensity(intensity, ndirs, nbands)
    n = size(intensity,2); # num cells
    omega = 2*pi/ndirs; # angle per direction*2  why *2?
    int_intensity = zeros(nbands, n);
    for i=1:n
        for j=1:nbands
            for k=1:ndirs
                int_intensity[j,i] += intensity[(j-1)*ndirs+k, i] * omega;
            end
        end
    end
    
    return int_intensity;
end

# Derivative of equilibrium intensity with respect to temperature
# input: frequency bands, bandwidth, temperature at each cell
# output: dI/dT at each cell
# NOTE: This produces an array of dIdT for every cell and band. Single cell/band version is below
function dIdT(freq, dw, temp; polarization="T")
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
    
    return didt;
end

# Derivative of equilibrium intensity with respect to temperature
# input: center frequency, bandwidth, temperature at one cell
# output: dI/dT for one cell and band
# NOTE: This is for one cell and band. Full cells/bands version above.
function dIdT_single(freq, dw, temp; polarization="T")
    if polarization=="T"
        vs = vs_TAS;
        c = c_TAS;
    else # L
        vs = vs_LAS;
        c = c_LAS;
    end
    
    tmp = 0.0;
    for gi=1:20 # gaussian quadrature
        fi = freq + dw/2 * g20xi[gi]; # frequency at gauss point
        # K2 = ((-vs + sqrt(vs*vs + 4*fi*c)) / (2*c))^2; # K^2
        tmpK = (-vs + sqrt(vs*vs + 4*fi*c)) / (2*c); # K updated to match ipcalc
        tmp2 = exp(hobol*fi/temp);
        # tmp += (tmp2 * fi * K2 / (tmp2 - 1)^2) * wi[gi];
        tmp += (tmp2 * (fi * tmpK)^2 / (tmp2 - 1)^2) * g20wi[gi]; # updated to match ipcalc
    end
    didt = tmp * dw * 3.2473482785757725e-48 / (temp * temp); # dirac * hobol / (8*pi^3) = 3.2473482785757725e-48
    
    return didt;
end

# Change temperature for one time step
# input: temperature array to update, previous step temperature, intensity from this and the previous time step, 
#        temperature from previous time step, frequency bands, bandwidth
# output: temperature for next step for each cell
function get_next_temp!(temp_next, temp_last, I_last, I_next, freq, dw; polarization="T")
    debug = false;
    n = length(temp_last); # number of cells
    m = length(freq); # number of bands
    
    # Io_last is equilibrium I from previous temp
    Io_last = Io.values;
    
    G_last = get_integrated_intensity(I_last, ndirs, nbands);
    G_next = get_integrated_intensity(I_next, ndirs, nbands);
    # didt = dIdT(freq, dw, temp_last, polarization=polarization); # Use single version in loop instead
    
    dt = Finch.time_stepper.dt;
    idt = 1/dt;
    
    if debug
        max_iters = 0;
        ave_iters = 0;
    end
    
    tol = 1e-10;
    maxiters = 50;
    for i=1:n # loop over cells
        # old values are not updated
        uold = 0.0; # These strange names are taken from the fortran code
        gnb = 0.0;  #
        for j=1:m # loop over bands
            uold += Io_last[j,i] * idt / group_v[j];
            gnb += G_last[j,i] * idt / group_v[j];
        end
        uold = 4*pi*uold;
        
        delta_T = 0;
        for iter=1:maxiters # iteratively refine delta_T
            unew = 0.0; # These strange names are taken from the fortran code
            gna = 0.0;  #
            uprime = 0.0;
            for j=1:m # loop over bands
                beta = 1 / get_time_scale(freq[j], temp_next[i], polarization=polarization);
                didt = dIdT_single(freq[j], dw, temp_last[i], polarization=polarization);
                
                unew += equilibrium_intensity(freq[j], dw, temp_next[i], polarization=polarization) * (beta + idt) / group_v[j];
                gna += G_next[j,i] * (beta + idt) / group_v[j];
                uprime += 4*pi * didt * (beta + idt) / group_v[j];
            end
            unew = 4*pi*unew;
            
            delta_T = (uold + gna - gnb - unew) / uprime;
            
            if debug println("cell "*string(i)*" ("*string(uold)*", "*string(unew)*", "*string(gna)*", "*string(gnb)*", "*string(uprime)*") : "*string(delta_T)) end
            
            temp_next[i] = temp_next[i] + delta_T;
            
            if abs(delta_T) < tol
                if debug
                    max_iters = max(max_iters, iter);
                    ave_iters += iter;
                end
                break;
            else
                if iter==maxiters
                    println("Temperature change didn't converge. Last delta_T = "*string(delta_T));
                end
            end
        end# iterative refinement
    end# cell loop
    if debug
        println("ave iterations: "*string(ave_iters/n)*"  max: "*string(max_iters))
    end
    
    return temp_next;
end

# Update temperature, equilibrium I, and time scale
#
function update_temperature(temp, I_last, I_next, freq, dw)
    
    temp_last = deepcopy(temp);
    
    get_next_temp!(temp, temp_last, I_last, I_next, freq, dw);
    
    equilibrium_intensity!(Io.values, freq, dw, temp);
    get_time_scale!(beta.values, freq, temp);
    
    # update I_last here as well
    for i=1:length(I_next)
        I_last[i] = I_next[i];
        if I_next[i] === NaN
            println("Error: NaN found in I");
            exit(0);
        end
    end
end

#############################################################################
## alternative models.

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