function equilibrium_intensity_noglobal(freq, dw, temp, polarization, g20xi, g20wi)::Float64
    vs::Float64 = 0.0;
    c::Float64 = 0.0;
    if polarization == 0.0
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

@callbackFunction(
    function isothermal_bdry_noglobal(
        I, vg, sx, sy, center_freq,
        polarizations, delta_freq, g20xi, g20wi,
        band, dir, normal, temp, ndir)
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2]

        if sdotn > 0 # outward
            interior_intensity::Float64 = I[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward gains from equilibrium
            center_f::Float64 = center_freq[band];
            polarization::Float64 = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            temp::Float64 = Float64(temp);
            iso_intensity::Float64 = equilibrium_intensity_noglobal(center_f, delta_f, temp, polarization, g20xi, g20wi);
            result = -vg[band] * iso_intensity * sdotn;
        end
        
        return result;
    end
)