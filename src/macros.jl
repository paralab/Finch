#=
Macros for the interface.
Many of these will eventually be removed.
Try to use the regular function interface in finch_interface.jl
=#
export @generateFor, @domain, @mesh, @solver, @stepper, @setSteps, @functionSpace, @trialFunction, @matrixFree,
        @testFunction, @nodes, @order, @boundary, @addBoundaryID, @referencePoint, @variable, @coefficient, 
        @parameter, @testSymbol, @initial, @preStepFunction, @postStepFunction, @callbackFunction,
        @timeInterval, @weakForm, @fluxAndSource, @LHS, @RHS, @exportCode, @importCode,
        @customOperator, @customOperatorFile,
        @outputMesh, @useLog, @finalize

###############################################################################
# These will be kept.

#= Assume something like this
myf = @callbackFunction(
    function f(a,b,c)
        body...
    end
)
=#
macro callbackFunction(f)
    fname = string(f.args[1].args[1]);
    fargs = [string(f.args[1].args[i]) for i=2:length(f.args[1].args)];
    fbody = string(f.args[2]);
    return esc(:(callbackFunction($f, name=$fname, args=$fargs, body=$fbody)));
end

macro preStepFunction(f)
    return esc(quote
        function PRE_STEP_FUNCTION()
            $f
        end
        preStepFunction(PRE_STEP_FUNCTION);
    end)
end

macro postStepFunction(f)
    return esc(quote
        function POST_STEP_FUNCTION()
            $f
        end
        postStepFunction(POST_STEP_FUNCTION);
    end)
end

###############################################################################
# These three could be kept for a minor convenience, but may be removed eventually.

macro variable(var) return esc(:(@variable($var, SCALAR, NODAL))); end
macro variable(var, type) return esc(:(@variable($var, $type, NODAL))); end
macro variable(var, type, location)
    varsym = string(var);
    return esc(quote
        $var = variable($varsym, $type, $location);
    end)
end

macro coefficient(c, val) return esc(:(@coefficient($c, SCALAR, $val))); end
macro coefficient(c, type, val)
    csym = string(c);
    return esc(quote
        coefficient($csym, $val, $type);
    end)
end

macro parameter(p,val) return esc(:(@parameter($p, SCALAR, $val))); end
macro parameter(p, type, val)
    psym = string(p);
    return esc(quote
        parameter($psym, $val, $type);
    end)
end

#################################################################################
# These below will eventually be removed

# Initializes code generation module
macro generateFor(lang) return esc(:(@generateFor($lang, Finch.project_name, ""))); end
macro generateFor(lang, filename) return esc(:(@generateFor($lang, $filename, ""))); end
macro generateFor(lang, filename, header)
    return esc(quote
        generateFor($lang, filename=$filename, header=$header);
    end)
end

# Optionally create a log file
macro useLog() return :(@uselog(Finch.project_name)) end
macro useLog(name) return esc(:(@useLog($name, Finch.output_dir))) end
macro useLog(name, dir)
    return esc(quote
        useLog($name, dir=$dir);
    end)
end

# Sets dimension, domain shape, discretization type
macro domain(dims) return :(@domain($dims, SQUARE, UNIFORM_GRID)) end
macro domain(dims, geometry, mesh)
    return esc(quote
        domain($dims, shape=$geometry, grid=$mesh);
    end)
end

# Set the solver type: CG, etc.
macro solver(type)
    return esc(quote
        solverType($type);
    end)
end

# Set basis function type and order
macro functionSpace(space) return esc(:(@functionSpace($space, 1, 1))) end
macro functionSpace(space, N) return esc(:(@functionSpace($space, $N, $N))) end
macro functionSpace(space, Nmin, Nmax)
    return esc(quote
        functionSpace(space=$space, orderMin=$Nmin, orderMax=$Nmax);
    end)
end

macro trialFunction(space, N) return esc(:(@trialFunction($space, $N, $N))) end
macro trialFunction(space, Nmin, Nmax)
    return esc(quote
        trialSpace(space=$space, orderMin=$Nmin, orderMax=Nmax);
    end)
end

macro testFunction(space, N) return esc(:(@testFunction($space, $N, $N))) end
macro testFunction(space, Nmin, Nmax)
    return esc(quote
        testSpace(space=$space, orderMin=$Nmin, orderMax=Nmax);
    end)
end

# Sets elemental node locations: Lobatto, etc.
macro nodes(nodetype)
    return esc(quote
        nodeType($nodetype);
    end)
end

# Sets time stepper type 
macro stepper(s) return esc(:(@stepper($s, 0))) end
macro stepper(s, cfl)
    return esc(quote
        timeStepper($s, cfl=$cfl);
    end)
end

# Sets time stepper dt and Nsteps to specified values
macro setSteps(dt, steps)
    return esc(quote
        setSteps($dt, $steps);
    end)
end

# Selects matrix free
macro matrixFree() return esc(:(@matrixFree(1000, 1e-6))) end
macro matrixFree(max, tol)
    return esc(quote
        matrixFree(maxiters=$max, tol=$tol);
    end)
end

# Adds a custom operator
macro customOperator(s, handle)
    symb = string(s);
    return esc(quote
        customOperator($symb, $handle);
    end)
end

# Adds a set of custom operators defined in a file
macro customOperatorFile(file)
    return esc(quote
        customOperatorFile($file);
    end)
end

# ============= end config , begin prob ==============

# This one-argument version should only be used for importing from a file
macro mesh(m)
    return esc(quote
        mesh($m);
    end)
end
# More than one argument is for generating a simple mesh.
macro mesh(m,N) return esc(quote @mesh($m,$N,1,[0,1]); end) end
macro mesh(m,N,bids) return esc(quote @mesh($m,$N,$bids,[0,1]); end) end
macro mesh(m, N, bids, interval)
    return esc(quote
        mesh($m, elsperdim=$N, bids=$bids, interval=$interval);
    end)
end

# Write the mesh to a MSH file
macro outputMesh(file) return esc(:(@outputMesh($file, MSH_V2))); end
macro outputMesh(file, format)
    return esc(quote
        exportMesh($file, $format);
    end)
end

macro addBoundaryID(bid, expression)
    return esc(quote
        trueOnBdry(x, y=0, z=0) = $expression; # points with x,y,z on this bdry segment evaluate true here
        addBoundaryID($bid, trueOnBdry);
    end)
end

macro boundary(var, bid, bc_type) return esc(:(@boundary($var,$bid,$bc_type, 0))); end
macro boundary(var, bid, bc_type, bc_exp)
    return esc(quote
        boundary($var, $bid, $bc_type, $bc_exp);
    end)
end

macro referencePoint(var, pos, val)
    return esc(quote
        referencePoint($var, $pos, $val);
    end)
end

macro timeInterval(t)
    return esc(quote
        timeInterval($t);
    end)
end

macro initial(var, ic)
    return esc(quote
        initial($var, $ic)
    end)
end

macro testSymbol(var) return esc(:(@testSymbol($var, SCALAR))); end
macro testSymbol(var, type)
    varsym = string(var);
    return esc(quote
        testSymbol($varsym, $type);
    end)
end

macro weakForm(var, ex)
    return esc(quote
        weakForm($var, $ex);
    end)
end

macro fluxAndSource(var, fex, sex)
    return esc(quote
        fluxAndSource($var, $fex, $sex);
    end)
end

# Import and export the code layer functions to allow manual changes
macro exportCode(file)
    return esc(quote
        exportCode($file);
    end)
end

macro importCode(file)
    return esc(quote
        importCode($file);
    end)
end

macro finalize()
    return esc(:(Finch.finalize()));
end
