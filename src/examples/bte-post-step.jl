function get_time_scale(freq, temp, polarization)
    if polarization == 0.0
        beta = btn * (temp*temp*temp*temp) * freq
        if freq >= wmax_half
            beta += btu * (freq*freq) / sinh(hobol * freq / temp);
        end
    else
        beta = bl * (temp*temp*temp) * (freq*freq);
    end

    return beta;
end

function dIdT_single(freq, dw, temp, polarization, g20xi, g20wi)
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization==0.0
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

# postStepFunction("update_temperature(temperature, temperatureLast, I, beta, center_freq, delta_freq, polarizations, G_last, G_next, vg, Io, g20xi, g20wi, ndirs, nbands, dt)")
@callbackFunction(
function update_temperature(temperature, temperatureLast, I, beta, center_freq, delta_freq, polarizations, G_last, G_next, vg, Io, g20xi, g20wi, ndirs, nbands, dt, cell)
    temperatureLast[cell] = temperature[cell]

    omega = 2*pi/ndirs * 2
    for band = 1:nbands
        G_next[band, cell] = 0.0
        for dir = 1:ndirs
            G_next[band, cell] += I[(band-1)*ndirs+dir, cell] * omega
        end
    end

    idt::Float64 = 1.0 / dt
    maxiters = 5

    uold = 0.0
    uchange = 0.0
    gnb = 0.0
    gna = 0.0
    for band = 1:nbands
        beta_f = get_time_scale(center_freq[band], temperature[cell], polarizations[band])

        uold += Io[band, cell] / vg[band]
        uchange += Io[band, cell] * beta_f / vg[band]
        gnb += G_last[band, cell] * (idt - beta_f) / vg[band]
        gna += G_next[band, cell] * idt / vg[band]

        G_last[band, cell] = G_next[band, cell]
    end
    uold *= 4*pi * idt
    uchange *= 4*pi

    unew = (uold + gna - gnb - uchange) / idt

    delta_T = 0
    for iter = 1:maxiters
         uchange = 0.0
         uprime = 0.0
         for band = 1:nbands
            didt = dIdT_single(center_freq[band], delta_freq[band], temperature[cell], polarizations[band], g20xi, g20wi)

            center_f::Float64 = center_freq[band]
            polarization::Float64 = polarizations[band]
            delta_f::Float64 = delta_freq[band]
            temp::Float64 = Float64(temperature[cell]);
            uchange += equilibrium_intensity_noglobal(center_f, delta_f, temp, polarization, g20xi, g20wi) / vg[band]
            uprime += didt / vg[band]
         end
         uchange *= 4*pi
         uprime *= 4*pi

         delta_u = unew - uchange
         delta_T = delta_u / uprime

         temperature[cell] = temperature[cell] + delta_T

         if iter == 1
            uold = delta_u
         end
    end

    for band = 1:nbands
        center_f::Float64 = center_freq[band]
        polarization::Float64 = polarizations[band]
        delta_f::Float64 = delta_freq[band]
        temp::Float64 = Float64(temperature[cell]);
        Io[band, cell] = equilibrium_intensity_noglobal(center_f, delta_f, temp, polarization, g20xi, g20wi)
        beta[band, cell] = get_time_scale(center_freq[band], temperature[cell], polarizations[band])
    end
end
)