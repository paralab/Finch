# The basic symbolic operators that are included automatically

function sym_dot_op(a,b)
    if (size(a) == size(b)) && (ndims(a) == 1)
        return [transpose(a)*b];
    else
        printerr("Usupported dimensions for: dot(a,b), sizes: a="*string(size(a))*", b="*string(size(b))*")");
    end
end

function sym_inner_op(a,b)
    if size(a) == size(b)
        c = 0;
        for i=1:length(a)
            c += a[i]*b[i];
        end
        return [c];
    else
        printerr("Unequal dimensions for: inner(a,b), sizes: a="*string(size(a))*", b="*string(size(b))*")");
    end
end

function sym_cross_op(a,b)
    if size(a) == size(b)
        if size(a) == (1,) # scalar
            # return [a[1]*b[1]]; # Should we allow this?
        elseif ndims(a) == 1 # vector
            if length(a) == 2 # 2D
                # return [a[1]*b[2] - a[2]*b[1]]; # Should we allow this?
            elseif length(a) == 3 # 3D
                return [a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]];
            end
        else
            printerr("Unsupported dimensions for: cross(a,b), sizes: a="*string(size(a))*", b="*string(size(b))*")");
        end
    else
        printerr("Unequal dimensions for: cross(a,b), sizes: a="*string(size(a))*", b="*string(size(b))*")");
    end
    printerr("Unexpected dimensions for: cross(a,b), sizes: a="*string(size(a))*", b="*string(size(b))*")");
    return [0];
end

function sym_transpose_op(a)
    if size(a) == (1,) # scalar
        return a;
    elseif ndims(a) == 1 # vector
        return reshape(a,1,length(a));
    else
        return permutedims(a);
    end
end

#########################################################################
# Special ops for other integrals

# Surface integral
function sym_surface_op(ex)
    newex = copy(ex);
    SURFACEINTEGRAL = symbols("SURFACEINTEGRAL");
    if typeof(ex) <: Array
        for i=1:length(ex)
            newex[i] = sym_surface_op(newex[i]);
        end
    elseif typeof(ex) == Basic
        return SURFACEINTEGRAL*ex;
    elseif typeof(ex) <: Number
        return SURFACEINTEGRAL*ex;
    end
    
    return newex;
end

# Boundary surface integral
function sym_boundary_op(ex)
    newex = copy(ex);
    BOUNDARYINTEGRAL = symbols("BOUNDARYINTEGRAL");
    if typeof(ex) <: Array
        for i=1:length(ex)
            newex[i] = sym_boundary_op(newex[i]);
        end
    elseif typeof(ex) == Basic
        return BOUNDARYINTEGRAL*ex;
    elseif typeof(ex) <: Number
        return BOUNDARYINTEGRAL*ex;
    end
    
    return newex;
end

# Dirichlet Boundary surface integral
function sym_dirichletBoundary_op(ex)
    newex = copy(ex);
    DIRICHLETINTEGRAL = symbols("DIRICHLETINTEGRAL");
    if typeof(ex) <: Array
        for i=1:length(ex)
            newex[i] = sym_dirichletBoundary_op(newex[i]);
        end
    elseif typeof(ex) == Basic
        return DIRICHLETINTEGRAL*ex;
    elseif typeof(ex) <: Number
        return DIRICHLETINTEGRAL*ex;
    end
    
    return newex;
end

# Neumann Boundary surface integral
function sym_neumannBoundary_op(ex)
    newex = copy(ex);
    NEUMANNINTEGRAL = symbols("NEUMANNINTEGRAL");
    if typeof(ex) <: Array
        for i=1:length(ex)
            newex[i] = sym_neumannBoundary_op(newex[i]);
        end
    elseif typeof(ex) == Basic
        return NEUMANNINTEGRAL*ex;
    elseif typeof(ex) <: Number
        return NEUMANNINTEGRAL*ex;
    end
    
    return newex;
end

# Special ops for DG

function sym_ave_op(var)
    #prefix = "DGAVERAGE_";
    side1 = "DGSIDE1_";
    side2 = "DGSIDE2_";
    if typeof(var) <: Array
        result = copy(var);
        for i=1:length(result)
            result[i] = sym_ave_op(var[i]);
        end
    elseif typeof(var) == Basic
        result = (symbols(side1*string(var)) + symbols(side2*string(var))) * Basic(0.5);
    elseif typeof(var) <: Number
        result = Basic(var);
    end
    return result;
end

function sym_jump_op(var)
    #prefix = "DGJUMP_";
    side1 = "DGSIDE1_";
    side2 = "DGSIDE2_";
    norm1 = sym_normal_op(1);
    norm2 = sym_normal_op(2);
    if typeof(var) <: Array
        # make into a vector according to normal
        if length(var) == 1
            result = sym_jump_op(var[1]);
        else
            # matrix rows for vec components, cols for normal components
            result = sym_jump_op(var[1]);
            for i=2:length(var)
                result = hcat(result, sym_jump_op(var[i]));
            end
            permutedims(result)
        end
        
    elseif typeof(var) == Basic
        # Note: norm is a vector
        result = (symbols(side1*string(var)) - symbols(side2*string(var))) .* norm1;
    elseif typeof(var) <: Number
        result = Basic(0) .* norm1;
    end
    
    return result;
end

# function sym_ave_normdotgrad_op(var)
#     prefix = "DGAVENORMDOTGRAD_";
#     if typeof(var) <: Array
#         result = copy(var);
#         for i=1:length(result)
#             result[i] = sym_ave_normdotgrad_op(var[i]);
#         end
#     elseif typeof(var) == Basic
#         result = symbols(prefix*string(var));
#     elseif typeof(var) <: Number
#         result = Basic(0);
#     end
#     return result;
# end

# function sym_jump_normdotgrad_op(var)
#     prefix = "DGJUMPNORMDOTGRAD_";
#     if typeof(var) <: Array
#         result = copy(var);
#         for i=1:length(result)
#             result[i] = sym_jump_normdotgrad_op(var[i]);
#         end
#     elseif typeof(var) == Basic
#         result = symbols(prefix*string(var));
#     elseif typeof(var) <: Number
#         result = Basic(0);
#     end
#     return result;
# end

function sym_normal_op()
    _FACENORMAL_1 = symbols("_FACENORMAL1_1");
    _FACENORMAL_2 = symbols("_FACENORMAL1_2");
    _FACENORMAL_3 = symbols("_FACENORMAL1_3");
    d = finch_state.config.dimension;
    if d==1
        return [_FACENORMAL_1];
    elseif d==2
        return [_FACENORMAL_1, _FACENORMAL_2];
    elseif d==3
        return [_FACENORMAL_1, _FACENORMAL_2, _FACENORMAL_3];
    end
end

function sym_normal_op(side)
    _FACENORMAL_1 = symbols("_FACENORMAL"*string(side)*"_1");
    _FACENORMAL_2 = symbols("_FACENORMAL"*string(side)*"_2");
    _FACENORMAL_3 = symbols("_FACENORMAL"*string(side)*"_3");
    d = finch_state.config.dimension;
    if d==1
        return [_FACENORMAL_1];
    elseif d==2
        return [_FACENORMAL_1, _FACENORMAL_2];
    elseif d==3
        return [_FACENORMAL_1, _FACENORMAL_2, _FACENORMAL_3];
    end
end

#############################
# Specific to shifted boundary method
function sym_trueNormal_op()
    _TRUENORMAL_1 = symbols("_TRUENORMAL_1");
    _TRUENORMAL_2 = symbols("_TRUENORMAL_2");
    _TRUENORMAL_3 = symbols("_TRUENORMAL_3");
    d = finch_state.config.dimension;
    if d==1
        return [_TRUENORMAL_1];
    elseif d==2
        return [_TRUENORMAL_1, _TRUENORMAL_2];
    elseif d==3
        return [_TRUENORMAL_1, _TRUENORMAL_2, _TRUENORMAL_3];
    end
end

function sym_distanceToBoundary_op()
    _DIST2BDRY_1 = symbols("_DIST2BDRY_1");
    _DIST2BDRY_2 = symbols("_DIST2BDRY_2");
    _DIST2BDRY_3 = symbols("_DIST2BDRY_3");
    d = finch_state.config.dimension;
    if d==1
        return [_DIST2BDRY_1];
    elseif d==2
        return [_DIST2BDRY_1, _DIST2BDRY_2];
    elseif d==3
        return [_DIST2BDRY_1, _DIST2BDRY_2, _DIST2BDRY_3];
    end
end

function sym_dirichletValue_op()
    return [symbols("BOUNDARYVALUE")];
end
function sym_neumannValue_op()
    return [symbols("BOUNDARYVALUE")];
end
function sym_elementDiameter_op()
    return [symbols("ELEMENTDIAMETER")];
end
function sym_elementVolume_op()
    return [symbols("ELEMENTVOLUME")];
end

#########################################################################
# derivative ops
#########################################################################

# Multiplies the expression by TIMEDERIV
# Dt(a*b+c)  ->  TIMEDERIV*(a*b+c)
function sym_Dt_op(ex)
    newex = copy(ex);
    TIMEDERIV = symbols("TIMEDERIV");
    if typeof(ex) <: Array
        for i=1:length(ex)
            newex[i] = sym_Dt_op(newex[i]);
        end
    elseif typeof(ex) == Basic
        return TIMEDERIV*ex;
    elseif typeof(ex) <: Number
        newex = 0;
    end
    
    return newex;
end

# Applies a derivative prefix. wrt is the axis index
# sym_deriv_string(u_12, 1) -> D1_u_12
function sym_deriv_string(var, wrt)
    # since numbers may have been converted to Basic[]
    if typeof(wrt) <: Number
        swrt = string(wrt);
    elseif typeof(wrt) <: Array
        swrt = string(wrt[1]);
    end
    
    prefix = "D"*swrt*"_";
    
    # var could be a single symbol like _u_1
    # or an expression like _u_1 + _v_1
    # To check, parse as Expr
    ex = Meta.parse(string(var));
    if typeof(ex) == Symbol
        result = symbols(prefix*string(var));
    elseif typeof(ex) <:Number
        result = 0;
    elseif typeof(ex) == Expr
        if ex.head === :call && (ex.args[1] in [:+ , :- , :.+ , :.-])
            for i=2:length(ex.args)
                ex.args[i] = sym_deriv_string(ex.args[i], wrt);
            end
            result = ex;
        elseif ex.head === :call && (ex.args[1] in [:* , :.*])
            newex = :(a+b);
            newex.args = [:+];
            for i=2:length(ex.args)
                subex = copy(ex);
                subex.args[i] = sym_deriv_string(ex.args[i], wrt);
                push!(newex.args, subex);
            end
            result = newex;
        else
            printerr("Unexpected expression in sym_deriv_string: "*string(var));
            result = 0;
        end
        result = Basic(result);
    else
        printerr("Unexpected expression in sym_deriv_string: "*string(var));
        result = 0;
    end
    
    return result;
end

function sym_deriv_op(u, wrt)
    if typeof(u) <: Array
        result = copy(u);
        for i=1:length(result)
            result[i] = sym_deriv_op(u[i], wrt);
        end
    elseif typeof(u) == Basic
        result = sym_deriv_string(u, wrt);
    elseif typeof(u) <: Number
        result = 0;
    end
    
    return result;
end
    
function sym_grad_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) <: Array
        d = finch_state.config.dimension;
        rank = 0;
        if ndims(u) == 1 && length(u) > 1
            rank = 1;
        elseif ndims(u) == 2 
            rank = 2;
        end
        
        if rank == 0
            # result is a vector
            for i=1:d
                push!(result, sym_deriv_string(u[1], i));
            end
        elseif rank == 1
            # result is a tensor
            for i=1:d
                for j=1:d
                    push!(result, sym_deriv_string(u[i], j));
                end
            end
            result = reshape(result, d,d);
        elseif rank == 2
            # not yet ready
            # printerr("unsupported operator, grad(tensor)");
            # return nothing;
            for i=1:d
                for j=1:d
                    for k=1:d
                        push!(result, sym_deriv_string(u[i+(j-1)*d], k));
                    end
                end
            end
            result = reshape(result, d,d,d);
        end
    elseif typeof(u) == Basic
        # result is a vector
        d = finch_state.config.dimension;
        for i=1:d
            push!(result, sym_deriv_string(u, i));
        end
    elseif typeof(u) <: Number
        return zeros(finch_state.config.dimension);
    end
    
    return result;
end

function sym_div_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) <: Array
        d = finch_state.config.dimension;
        rank = 0;
        if ndims(u) == 1 && length(u) > 1
            rank = 1;
        elseif ndims(u) == 2 
            rank = 2;
        end
        
        if rank == 0
            # Not allowed
            printerr("unsupported operator, div(scalar)");
            return nothing;
        elseif rank == 1
            # result is a scalar
            if d==1
                result = [sym_deriv_string(u[1], 1)];
            else
                ex = :(a+b);
                ex.args = [:+];
                for i=1:d
                    push!(ex.args, sym_deriv_string(u[i], i))
                end
                result = [Basic(ex)];
            end
        elseif rank == 2
            # not yet ready
            printerr("unsupported operator, div(tensor)");
            return nothing;
        end
    elseif typeof(u) <: Number
        # Not allowed
        printerr("unsupported operator, div(number)");
        return nothing;
    end
    
    return result;
end

function sym_curl_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) <: Array
        d = finch_state.config.dimension;
        rank = 0;
        if ndims(u) == 1 && sz[1] > 1
            rank = 1;
        elseif ndims(u) == 2 
            rank = 2;
        end
        
        if rank == 0
            # Not allowed
            printerr("unsupported operator, curl(scalar)");
            return nothing;
        elseif rank == 1
            # result is a vector
            if d==1
                result = [sym_deriv_string(u[1], 1)];
            else
                #TODO
                printerr("curl not ready");
                return nothing;
            end
        elseif rank == 2
            # not yet ready
            printerr("unsupported operator, curl(tensor)");
            return nothing;
        end
    elseif typeof(u) <: Number
        # Not allowed
        printerr("unsupported operator, curl(number)");
        return nothing;
    end
    
    return result;
end

function sym_laplacian_op(u)
    # simply use the above ops
    return sym_div_op(sym_grad_op(u));
end

#################################################################################
# The finite volume operators.
function sym_left_op(f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_left_op(f[i]);
        end
        
    elseif typeof(f) == Basic
        # Apply a CELL1 tag to signal to the code generator which side of the face the value is from.
        result = apply_flag_to_all_symbols("CELL1", f);
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_right_op(f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_right_op(f[i]);
        end
        
    elseif typeof(f) == Basic
        # Apply a CELL2 tag to signal to the code generator which side of the face the value is from.
        result = apply_flag_to_all_symbols("CELL2", f);
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_central_op(f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_right_op(f[i]);
        end
        
    elseif typeof(f) == Basic
        # Apply a CELL2 tag to signal to the code generator which side of the face the value is from.
        result = apply_flag_to_all_symbols("CENTRAL", f);
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_neighborhood_op(f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_neighborhood_op(f[i]);
        end
        
    elseif typeof(f) == Basic
        # Apply a CELL2 tag to signal to the code generator which side of the face the value is from.
        result = apply_flag_to_all_symbols("NEIGHBORHOOD", f);
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_upwind_op(v, f)
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_upwind_op(v, f[i]);
        end
        
    elseif typeof(f) == Basic
        @funs(conditional)
        @funs(isgreaterthan)
        # Apply a CELLn tag to signal to the code generator which side of the face the value is from.
        side1 = apply_flag_to_all_symbols("CELL1", f);
        side2 = apply_flag_to_all_symbols("CELL2", f);
        # The formula for an adjustable first order upwind scheme.
        # F(u) = 0.5*(side1+side2)*(v.normal) + 0.5*(side1-side2)*abs(v.normal)*(1-alpha)
        result = conditional(isgreaterthan(sym_dot_op(v, sym_normal_op())[1], Basic(0)), 
                    (sym_dot_op(v, sym_normal_op()) .* side1)[1],
                    # else
                    (sym_dot_op(v, sym_normal_op()) .* side2)[1]) # Note that the arguments to conditional and isgreaterthan are not in arrays.
        
    elseif typeof(f) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(f);
    end
    return result;
end

function sym_upwindA_op(v, f, alpha=Basic(0))
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(f) <: Array
        result = copy(f);
        for i=1:length(result)
            result[i] = sym_upwindA_op(v, f[i], alpha);
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
        @funs(conditional)
        @funs(isgreaterthan)
        @funs(symbolmin)
        @funs(symbolmax)
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

###################################################################
# These are basic math ops like sin and exp
function sym_exp_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_exp_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(exp)
        result = exp(e);
    else #if typeof(e) <: Number
        result = Basic(exp(e));
    end
    return result;
end

function sym_sin_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_sin_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(sin)
        result = sin(e);
    else #if typeof(e) <: Number
        result = Basic(sin(e));
    end
    return result;
end

function sym_cos_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_cos_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(cos)
        result = cos(e);
    else #if typeof(e) <: Number
        result = Basic(cos(e));
    end
    return result;
end

function sym_tan_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_tan_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(tan)
        result = tan(e);
    else #if typeof(e) <: Number
        result = Basic(tan(e));
    end
    return result;
end

function sym_abs_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_abs_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(abs)
        result = abs(e);
    else #if typeof(e) <: Number
        result = Basic(abs(e));
    end
    return result;
end

function sym_sinh_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_sinh_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(sinh)
        result = sinh(e);
    else #if typeof(e) <: Number
        result = Basic(sinh(e));
    end
    return result;
end

function sym_cosh_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_cosh_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(cosh)
        result = cosh(e);
    else #if typeof(e) <: Number
        result = Basic(cosh(e));
    end
    return result;
end

function sym_tanh_op(e)
    if typeof(e) <: Array
        result = copy(e);
        for i=1:length(result) result[i] = sym_tanh_op(e[i]); end
    elseif typeof(e) == Basic
        @funs(tanh)
        result = tanh(e);
    else #if typeof(e) <: Number
        result = Basic(tanh(e));
    end
    return result;
end


# Load them into the global arrays
function load_basic_ops()
    op_names = [:dot, :inner, :cross, :transpose, :surface, :boundary, :dirichletBoundary, :neumannBoundary, :ave, :jump, 
                :normal, :trueNormal, :distanceToBoundary, :dirichletValue, :neumannValue, :elementDiameter, :elementVolume,
                :Dt, :deriv, :grad, :div, :curl, :laplacian,
                :left, :right, :central, :neighborhood, :upwind, :upwindA, :burgerGodunov,
                :exp, :sin, :cos, :tan, :abs, :sinh, :cosh, :tanh];
    _handles = [sym_dot_op, sym_inner_op, sym_cross_op, sym_transpose_op, sym_surface_op, sym_boundary_op, sym_dirichletBoundary_op,
                sym_neumannBoundary_op, sym_ave_op, sym_jump_op, 
                sym_normal_op, sym_trueNormal_op, sym_distanceToBoundary_op, sym_dirichletValue_op, sym_neumannValue_op,
                sym_elementDiameter_op, sym_elementVolume_op,
                sym_Dt_op, sym_deriv_op, sym_grad_op, 
                sym_div_op, sym_curl_op, sym_laplacian_op,
                sym_left_op, sym_right_op, sym_central_op, sym_neighborhood_op, sym_upwind_op, sym_upwindA_op, sym_burgerGodunov_op,
                sym_exp_op, sym_sin_op, sym_cos_op, sym_tan_op, sym_abs_op, sym_sinh_op, sym_cosh_op, sym_tanh_op];
    ops = Vector{SymOperator}(undef,0);
    for i=1:length(op_names)
        push!(ops, SymOperator(op_names[i], _handles[i]));
    end
    
    return ops;
end
