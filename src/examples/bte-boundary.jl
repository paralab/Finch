@callbackFunction(
    function adiabatic_bdry(intensity, vg, sx, sy, dir, band, normal)
        # debug = (fid == 1 && t > 0.29);
        result = intensity[dir + (band-1)*ndirs];
        
        specularity = 0.5;
        #specular_tol = 1e-5; # tolerance for selecting specular reflection direction
        reflect_tol = 1e-6; # tolerance for ignoring reflection (parallel to bdry -> no reflection)
        s = [sx[dir], sy[dir]];
        s_dot_n = s[1] * normal[1] + s[2] * normal[2];
        # if debug println("dir="*string(dir)*", normal="*string(normal)*", s_dot_n="*string(s_dot_n)*", val="*string(result)) end
        
        # Does the direction point in or out of the domain at this boundary?
        if s_dot_n < -1e-6 # into the interior
            sspec = s .- 2*s_dot_n .* normal; # specular reflection vector
            close_to_specular = 0.5; 
            angle = 2*pi/ndirs;
            
            # For every direction that points out of the boundary, add some reflected intensity to this.
            specular_part = 0;
            diffuse_part = 0;
            ndiffuse = 0;
            for di=1:ndirs
                si = [sx[di], sy[di]];
                si_dot_n = si[1]*normal[1] + si[2]*normal[2];
                if si_dot_n > 1e-6 # out of bdry and not too close to parallel
                    # Specular (takes dos portion of reflection directly)
                    #how_close = abs(si[1] - sspec[1])+abs(si[2] - sspec[2]);
                    #if how_close < close_to_specular
                    if reflect[dir] == di
                        specular_part = specularity * intensity[di + (band-1)*ndirs] * si_dot_n;
                        close_to_specular = how_close;
                        # if debug println("specular got "*string(specular_part)*" from "*string(di)) end
                    end
                    
                    # Diffuse (takes its share of the (1-dos) portion of reflection)
                    diffuse_part += (1-specularity) * intensity[di + (band-1)*ndirs] * si_dot_n;
                elseif si_dot_n < -1e-6
                    ndiffuse += 1; # count directions taking a share
                end
            end
            # Adjust diffuse part share depending on number of directions taking a share
            diffuse_part /= ndiffuse;
            # if debug println("diffuse got "*string(diffuse_part)*" from "*string(ndiffuse)*" dirs") end
            result = vg[band] * (specular_part + diffuse_part);
        elseif s_dot_n > 1e-6 # out
            result = -vg[band] * intensity[dir + (band-1)*ndirs] * s_dot_n;
            # if debug println("lost "*string(intensity[dir + (band-1)*ndirs] * s_dot_n)) end
        else
            result = 0;
        end
        
        return result;
    end
)
@callbackFunction(
    function isothermal_bdry(intensity, vg, sx, sy, band, dir, normal, temp)
        interior_intensity = intensity[dir + (band-1)*ndirs];
        iso_intensity = equilibrium_intensity(center_freq[band], delta_freq, temp);
        sdotn = sx[dir]*normal[1] + sy[dir]*normal[2];
        if sdotn > 0
            result = -vg[band] * interior_intensity * sdotn;
        else
            result = -vg[band] * iso_intensity * sdotn;
        end
        return result;
    end
)
@callbackFunction(
    function symmetric_bdry(intensity, vg, sx, sy, band, dir, normal)
        # This is essentially the same as pure specular reflection?
        interior_intensity = intensity[dir + (band-1)*ndirs];
        sym_intensity = intensity[reflect[dir] + (band-1)*ndirs];
        sdotn = sx[dir]*normal[1] + sy[dir]*normal[2];
        if sdotn > 0
            result = -vg[band] * interior_intensity * sdotn;
        else
            result = -vg[band] * sym_intensity * sdotn;
        end
        return result;
    end
)