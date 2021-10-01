# The finite volume operators that are included automatically if FV is used.

function sym_central_op(f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_central_op(f[i]);
        end
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
        
    elseif typeof(f) == Basic
        # Apply a CELLn tag to signal to the code generator which side of the face the value is from.
        side1 = apply_flag_to_all_symbols("CELL1", f);
        side2 = apply_flag_to_all_symbols("CELL2", f);
        result = Basic(0.5) .* (side1 .+ side2);
    end
    return result;
end

function sym_upwind_op(v, f, alpha=Basic(0))
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_upwind_op(v, f[i], alpha);
        end
        
    elseif typeof(f) == Basic
        # Apply a CELLn tag to signal to the code generator which side of the face the value is from.
        side1 = apply_flag_to_all_symbols("CELL1", f);
        side2 = apply_flag_to_all_symbols("CELL2", f);
        # The formula for an adjustable first order upwind scheme.
        # F(u) = 0.5*(side1+side2)*(v.normal) + 0.5*(side1-side2)*abs(v.normal)*(1-alpha)
        result = Basic(0.5) .* ( (sym_dot_op(v, sym_normal_op())) .* (side1 + side2) .+ abs.(sym_dot_op(v, sym_normal_op())[1]) .* (Basic(1.0)-alpha) .* (side1 - side2) );
        result = result[1]; # Since dot and normal returned arrays, result was placed in an array, but this is done above redundantly.
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_burgerGodunov_op(u, f)
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_burgerGodunov_op(u[i], f[i]);
        end
    elseif typeof(f) == Basic
        uside1 = apply_flag_to_all_symbols("CELL1", u);
        uside2 = apply_flag_to_all_symbols("CELL2", u);
        fside1 = apply_flag_to_all_symbols("CELL1", f);
        fside2 = apply_flag_to_all_symbols("CELL2", f);
        result = conditional(isgreaterthan(sym_normal_op()[1], Basic(0)), # if normal[1] > 0
                    conditional(isgreaterthan(uside1, uside2), # if u1 < u2
                        symbolmax(fside1, fside2), # max(f1,f2)
                        # else
                        conditional(isgreaterthan(uside1*uside2, Basic(0)), # if u1*u2 > 0 (same sign)
                            symbolmin(fside1, fside2), # min(f1,f2)
                            # else
                            Basic(0))), # 0
                    # elseif normal < 0 (just reverse 1 and 2)
                    conditional(isgreaterthan(uside2, uside1), # if u2 < u1
                        symbolmax(fside2, fside1), # max(f1,f2)
                        # else
                        conditional(isgreaterthan(uside2*uside1, Basic(0)), # if u1*u2 > 0 (same sign)
                            symbolmin(fside2, fside1), # min(f1,f2)
                            # else
                            Basic(0)))) # 0
    elseif typeof(f) <: Number
        result = Basic(f);
    end
    return result;
end

_names = [:upwind, :burgerGodunov, :central];
_handles = [sym_upwind_op, sym_burgerGodunov_op, sym_central_op];
