# Isothermal boundary assumes equilibrium intensity outside
@callbackFunction(
    function isothermal_bdry_axi(intensity, vg::Vector, sx::Vector, sy::Vector, isx::Vector, isy::Vector, 
                            band::Int, dir::Int, omega::Vector, normal::Vector{Float64}, temp)
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        insdotn::Float64 = isx[dir]*normal[1] + isy[dir]*normal[2];
        invomega::Float64 = 1.0/omega[dir];
        
        if sdotn > 0 # outward
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * insdotn;
            
        else # inward gains from equilibrium
            center_f::Float64 = center_freq[band];
            polarization::String = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            temp::Float64 = Float64(temp);
            iso_intensity::Float64 = equilibrium_intensity(center_f, delta_f, temp, polarization);
            result = -vg[band] * iso_intensity * insdotn;
        end
        
        # if band == 1 && dir == 1
        #     println("isothermal bdry sdotn $sdotn result $result")
        # end
        
        return result*invomega;
    end
)

# Symmetric boundary assumes mirror image of boundary cells
@callbackFunction(
    function symmetric_bdry_axi(intensity, vg::Vector, sx::Vector, sy::Vector, sz::Vector, isx::Vector, isy::Vector, 
                            band::Int, dir::Int, omega::Vector, normal::Vector{Float64})
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        insdotn::Float64 = isx[dir]*normal[1] + isy[dir]*normal[2];
        invomega::Float64 = 1.0/omega[dir];
        
        if sdotn > 0 # outward
            # use interior intensity
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * insdotn;
            
        else # inward
            # Find the reflection direction
            # Reflected vector is S - 2*Sdotn*n
            reflect_x = sx[dir] - 2*sdotn * normal[1];
            reflect_y = sy[dir] - 2*sdotn * normal[2];
            reflect_z = pi - sz[dir];
            closest = 1;
            difference = 0.0
            
            for i=1:ndir
                tmp = sx[i]*reflect_x + sy[i]*reflect_y + sz[i]*reflect_z;
                if tmp > difference
                    closest = i;
                    difference = tmp;
                end
            end
            
            sym_intensity::Float64 = intensity[closest + (band-1)*ndir];
            result = -vg[band] * sym_intensity * insdotn;
        end
        return result*invomega;
    end
)

# ALSI boundary is essentially the same as isothermal, but temperature is from the transducer.
@callbackFunction(
    function alsi_bdry_axi(intensity, vg::Vector, sx::Vector, sy::Vector, isx::Vector, isy::Vector, 
                            band::Int, dir::Int, omega::Vector, normal::Vector{Float64}, fid::Int)
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        insdotn::Float64 = isx[dir]*normal[1] + isy[dir]*normal[2];
        invomega::Float64 = 1.0/omega[dir];
        
        if sdotn > 0 # outward
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * insdotn;
            
        else # inward gains from equilibrium at transducer temp
            center_f::Float64 = center_freq[band];
            polarization::String = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            # find the index in alsi faces
            alsi::ALSIInfo = alsi_info;
            nfaces::Int = alsi.nfaces;
            alsifaces::Vector{Int} = alsi.faces;
            alsiindex = 0;
            for i=1:nfaces
                if alsifaces[i] == fid
                    alsiindex = i;
                    break;
                end
            end
            temp::Float64 = alsi.temp[alsiindex];
            iso_intensity::Float64 = equilibrium_intensity(center_f, delta_f, temp, polarization);
            result = -vg[band] * iso_intensity * insdotn;
        end
        
        return result*invomega;
    end
)

#############################################################################################################
# Below are regular, non-axisymetric BCs

# This is a 2D reflective (adiabatic) boundary
@callbackFunction(
    function adiabatic_bdry(intensity, vg::Vector, sx::Vector, sy::Vector, dir::Int, band::Int, normal::Vector{Float64})
        # debug = (fid == 1 && t > 0.29);
        result::Float64 = intensity[dir + (band-1)*ndirs];
        
        specularity = 0.5;
        reflect_tol = 1e-6; # tolerance for ignoring reflection (parallel to bdry -> no reflection)
        
        sdotn = sx[dir] * normal[1] + sy[dir] * normal[2];
        
        # Does the direction point in or out of the domain at this boundary?
        if sdotn < reflect_tol # into the interior
            # Find the closest direction vector to reflection
            reflect_x = sx[dir] - 2*sdotn * normal[1];
            reflect_y = sy[dir] - 2*sdotn * normal[2];
            reflect_dir = 1;
            difference = 0.0
            for i=1:ndir
                tmp = sx[i]*reflect_x + sy[i]*reflect_y;
                if tmp > difference
                    reflect_dir = i;
                    difference = tmp;
                end
            end
            
            # For every direction that points out of the boundary, add some reflected intensity to this.
            specular_part = 0;
            diffuse_part = 0;
            ndiffuse = 0;
            for di=1:ndirs
                si = [sx[di], sy[di]];
                si_dot_n = si[1]*normal[1] + si[2]*normal[2];
                if si_dot_n > reflect_tol # out of bdry and not too close to parallel
                    # Specular (takes dos portion of reflection directly)
                    if reflect_dir == di
                        specular_part = specularity * intensity[di + (band-1)*ndirs] * si_dot_n;
                    end
                    
                    # Diffuse (takes its share of the (1-dos) portion of reflection)
                    diffuse_part += (1-specularity) * intensity[di + (band-1)*ndirs] * si_dot_n;
                    
                elseif si_dot_n < -reflect_tol
                    ndiffuse += 1; # count directions taking a share
                end
            end
            # Adjust diffuse part share depending on number of directions taking a share
            diffuse_part /= ndiffuse;
            result = vg[band] * (specular_part + diffuse_part);
            
        elseif sdotn > 1e-6 # out of boundary
            result = -vg[band] * intensity[dir + (band-1)*ndirs] * sdotn;
        else
            result = 0;
        end
        
        return result;
    end
)

# Isothermal boundary assumes equilibrium intensity outside
@callbackFunction(
    function isothermal_bdry(intensity, vg::Vector, sx::Vector, sy::Vector, 
                            band::Int, dir::Int, normal::Vector{Float64}, temp)
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        
        if sdotn > 0 # outward
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward gains from equilibrium
            center_f::Float64 = center_freq[band];
            polarization::String = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            temp::Float64 = Float64(temp);
            iso_intensity::Float64 = equilibrium_intensity(center_f, delta_f, temp, polarization);
            result = -vg[band] * iso_intensity * sdotn;
        end
        
        return result;
    end
)

# Symmetric boundary assumes mirror image of boundary cells
@callbackFunction(
    function symmetric_bdry(intensity, vg::Vector, sx::Vector, sy::Vector, 
                            band::Int, dir::Int, normal::Vector{Float64})
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        if sdotn > 0 # outward
            # use interior intensity
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward
            # Find the reflection direction
            # Reflected vector is S - 2*Sdotn*n
            reflect_x = sx[dir] - 2*sdotn * normal[1];
            reflect_y = sy[dir] - 2*sdotn * normal[2];
            closest = 1;
            difference = 0.0
            
            for i=1:ndir
                tmp = sx[i]*reflect_x + sy[i]*reflect_y;
                if tmp > difference
                    closest = i;
                    difference = tmp;
                end
            end
            
            sym_intensity::Float64 = intensity[closest + (band-1)*ndir];
            result = -vg[band] * sym_intensity * sdotn;
        end
        return result;
    end
)

# ALSI boundary is essentially the same as isothermal, but temperature is from the transducer.
@callbackFunction(
    function alsi_bdry(intensity, vg::Vector, sx::Vector, sy::Vector, 
                            band::Int, dir::Int, normal::Vector{Float64})
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2];
        
        if sdotn > 0 # outward
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward gains from equilibrium at transducer temp
            center_f::Float64 = center_freq[band];
            polarization::String = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            temp::Float64 = 300.0; ############# TODO
            iso_intensity::Float64 = equilibrium_intensity(center_f, delta_f, temp, polarization);
            result = -vg[band] * iso_intensity * sdotn;
        end
        
        return result;
    end
)

# 3D versions of isothermal and symmetric boundaries
@callbackFunction(
    function isothermal_bdry_3d(intensity, vg::Vector, sx::Vector, sy::Vector, sz::Vector,
                                band::Int, dir::Int, normal::Vector{Float64}, temp)
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2] + sz[dir]*normal[3];
        
        if sdotn >= 0
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward
            center_f::Float64 = center_freq[band];
            polarization::String = polarizations[band];
            delta_f::Float64 = delta_freq[band];
            temp::Float64 = Float64(temp);
            iso_intensity::Float64 = equilibrium_intensity(center_f, delta_f, temp, polarization);
            result = -vg[band] * iso_intensity * sdotn;
        end
        
        return result;
    end
)
@callbackFunction(
    function symmetric_bdry_3d(intensity, vg::Vector, sx::Vector, sy::Vector, sz::Vector,
                                band::Int, dir::Int, normal::Vector{Float64})
        #
        ndir::Int = ndirs;
        sdotn::Float64 = sx[dir]*normal[1] + sy[dir]*normal[2] + sz[dir]*normal[3];
        if sdotn > 0 # outward
            # use interior intensity
            interior_intensity::Float64 = intensity[dir + (band-1)*ndir];
            result = -vg[band] * interior_intensity * sdotn;
            
        else # inward
            # Find the reflection direction
            # Reflected vector is S - 2*Sdotn*n
            reflect_x = sx[dir] - 2*sdotn * normal[1];
            reflect_y = sy[dir] - 2*sdotn * normal[2];
            reflect_z = sz[dir] - 2*sdotn * normal[3];
            reflect_dir = 1;
            difference = 0.0
            for i=1:ndir
                tmp = sx[i]*reflect_x + sy[i]*reflect_y + sz[i]*reflect_z;
                if tmp > difference
                    reflect_dir = i;
                    difference = tmp;
                end
            end
            
            sym_intensity::Float64 = intensity[reflect_dir + (band-1)*ndir];
            result = -vg[band] * sym_intensity * sdotn;
        end
        
        return result;
    end
)