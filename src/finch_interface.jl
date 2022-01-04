#=
This file contains all of the common interface functions.
Many of them simply call corresponding functions in jl.
=#
export generateFor, useLog, domain, solverType, functionSpace, trialSpace, testSpace, finiteVolumeOrder,
        nodeType, timeStepper, setSteps, matrixFree, customOperator, customOperatorFile,
        mesh, exportMesh, variable, coefficient, parameter, testSymbol, index, boundary, addBoundaryID,
        referencePoint, timeInterval, initial, preStepFunction, postStepFunction, callbackFunction,
        variableTransform, transformVariable,
        weakForm, fluxAndSource, flux, source, assemblyLoops,
        exportCode, importCode, printLatex,
        solve, cachesimSolve, finalize_finch, cachesim, output_values,
        morton_nodes, hilbert_nodes, tiled_nodes, morton_elements, hilbert_elements, 
        tiled_elements, ef_nodes, random_nodes, random_elements

# Begin configuration setting functions

function generateFor(lang; filename=project_name, header="", params=nothing)
    outputDirPath = pwd()*"/"*filename;
    if !isdir(outputDirPath)
        mkdir(outputDirPath);
    end
    framew = 0;
    if !in(lang, [CPP,MATLAB,DENDRO,HOMG])
        # lang should be a filename for a custom target
        # This file must include these three functions:
        # 1. get_external_language_elements() - file extensions, comment chars etc.
        # 2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
        # 3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code
        include(lang);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    else
        # Use an included target
        target_dir = @__DIR__
        if lang == DENDRO
            framew = DENDRO;
            lang = CPP;
            target_file = "/targets/target_dendro_cg.jl";
        elseif lang == MATLAB
            framew = MATLAB;
            lang = MATLAB;
            target_file = "/targets/target_matlab_cg.jl";
        else
            framew = 0;
            target_file = "/targets/target_matlab_cg.jl";
        end
        include(target_dir * target_file);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    end
    if !(params === nothing)
        set_codegen_parameters(params);
    end
end

function useLog(name=project_name; dir=output_dir, level=2)
    init_log(name, dir, level);
end

function domain(dims; shape=SQUARE, grid=UNIFORM_GRID)
    config.dimension = dims;
    config.geometry = shape;
    config.mesh_type = grid;
end

function solverType(type)
    set_solver(type);
end

function functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    config.trial_function = space;
    config.test_function = space;
    if orderMax > orderMin && orderMin >= 0
        config.p_adaptive = true;
        config.basis_order_min = orderMin;
        config.basis_order_max = orderMax;
    else
        config.basis_order_min = max(order, orderMin);
        config.basis_order_max = max(order, orderMin);
    end
end

function trialSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

function testSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

function nodeType(type)
    config.elemental_nodes = type;
end

function timeStepper(type; cfl=0)
    set_stepper(type, cfl);
end

function setSteps(dt, steps)
    set_specified_steps(dt, steps);
end

function matrixFree(shallwe=true; maxiters=100, tol=1e-6)
    config.linalg_matrixfree = shallwe;
    config.linalg_matfree_max = maxiters;
    config.linalg_matfree_tol = tol;
end

function customOperator(name, handle)
    s = Symbol(name);
    add_custom_op(s, handle);
end

function customOperatorFile(filename)
    log_entry("Adding custom operators from file: "*string(filename), 2);
    add_custom_op_file(filename);
end

# End configuration functions, begin problem definition functions

function mesh(msh; elsperdim=5, bids=1, interval=[0,1])
    if msh == LINEMESH
        log_entry("Building simple line mesh with nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_line_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    elseif msh == QUADMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple quad mesh with nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_quad_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    elseif msh == HEXMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_hex_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    else # msh should be a mesh file name
        # open the file and read the mesh data
        mfile = open(msh, "r");
        log_entry("Reading mesh file: "*msh);
        meshtime = @elapsed(mshdat=read_mesh(mfile));
        log_entry("Mesh reading took "*string(meshtime)*" seconds");
        close(mfile);
        # Assume an irregular mesh
        config.geometry = IRREGULAR;
        config.mesh_type = UNSTRUCTURED;
    end
    
    # Set this in Finch and build the corresponding Grid
    add_mesh(mshdat);
    
    # If bids>1 were specified for built meshes, add them here
    if bids > 1
        if msh == LINEMESH
            # already done in mesh_data
            # if bn == 2
            #     add_boundary_ID_to_grid(2, x -> (x >= interval[2]), grid_data);
            # end
            
        elseif msh == QUADMESH
            if bids == 2
                add_boundary_ID_to_grid(2, (x,y) -> (y <= interval[3]) || (y >= interval[4]), grid_data);
            elseif bids == 3
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> ((y <= interval[3]) || (y >= interval[4])) && (x > interval[1] && x < interval[2]), grid_data);
            elseif bids == 4
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> (y <= interval[3] && (x > interval[1] && x < interval[2])), grid_data);
                add_boundary_ID_to_grid(4, (x,y) -> (y >= interval[4] && (x > interval[1] && x < interval[2])), grid_data);
            end
            
        elseif msh == HEXMESH
            tiny = 1e-13;
            if bids == 6
                # bids = [1,2,3,4,5,6]; # all separate
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny), grid_data);
                add_boundary_ID_to_grid(6, (x,y,z) -> (z >= interval[6]-tiny), grid_data);
            elseif bids == 5
                # bids = [1,2,3,4,5]; # combine z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            elseif bids == 4
                # bids = [1,2,3,4]; # combine y and z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> ((y <= interval[3]+tiny) || (y >= interval[4]-tiny)), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> ((z <= interval[5]+tiny) || (z >= interval[6]-tiny)), grid_data);
            elseif bids == 3
                # bids = [1,2,3]; # combine x,y,z
                add_boundary_ID_to_grid(2, (x,y,z) -> (y <= interval[3]+tiny) || (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            elseif bids == 2
                # bids = [1,2]; # x=0, other
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny) || (y <= interval[3]+tiny) || (y >= interval[4]-tiny) || (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            end
            
        end
    end
end

function exportMesh(filename, format=MSH_V2)
    # open the file to write to
    mfile = open(filename, "w");
    log_entry("Writing mesh file: "*filename);
    output_mesh(mfile, format);
    close(mfile);
end

function finiteVolumeOrder(order)
    if config.dimension > 2
        printerr("Sorry, higher order FV is not ready for 3D (TODO: build parent/child grid)\n Continuing with first order.");
        return;
    end
    (parent, child) = divide_parent_grid(grid_data, order);
    set_parent_and_child(parent, child, order);
end

function variable(name, type=SCALAR, location=NODAL; index=nothing)
    varind = var_count + 1;
    varsym = Symbol(name);
    # Just make an empty variable with the info known so far.
    var = Variable(varsym, [], varind, type, location, [], index, 0, [], false);
    add_variable(var);
    return var;
end

function coefficient(name, val, type=SCALAR, location=NODAL; element_array=false)
    csym = Symbol(name);
    nfuns = makeFunctions(val); # if val is constant, nfuns will be 0
    return add_coefficient(csym, type, location, val, nfuns, element_array);
end

function parameter(name, val, type=SCALAR)
    if length(parameters) == 0
        coefficient("parameterCoefficientForx", "x")
        coefficient("parameterCoefficientFory", "y")
        coefficient("parameterCoefficientForz", "z")
        coefficient("parameterCoefficientFort", "t")
    end
    if typeof(val) <: Number
        newval = [val];
    elseif typeof(val) == String
        # newval will be an array of expressions
        # search for x,y,z,t symbols and replace with special coefficients like parameterCoefficientForx
        newval = [swap_parameter_xyzt(Meta.parse(val))];
        
    elseif typeof(val) <: Array
        newval = Array{Expr,1}(undef,length(val));
        for i=1:length(val)
            newval[i] = swap_parameter_xyzt(Meta.parse(val[i]));
        end
        newval = reshape(newval,size(val));
    else
        println("Error: use strings to define parameters");
        newval = 0;
    end
    
    return add_parameter(Symbol(name), type, newval);
end

function testSymbol(symb, type=SCALAR)
    add_test_function(Symbol(symb), type);
end

function index(name; range=[1])
    if length(range) == 2
        range = Array(range[1]:range[2]);
    end
    idx = Indexer(Symbol(name), range, range[1])
    add_indexer(idx);
    return idx;
end

function variableTransform(var1, var2, func)
    # Make sure things are valid
    if typeof(var1) <: Array
        if !(typeof(var2) <: Array)
            return add_variable_transform(var1, [var2], func);
        else
            return add_variable_transform(var1, var2, func);
        end
    else
        if typeof(var2) <: Array
            return add_variable_transform([var1], var2, func);
        else
            return add_variable_transform(var1, var2, func);
        end
    end
end

function transformVariable(xform::VariableTransform)
    transform_variable_values(xform);
end
function transformVariable(var1, var2)
    # Look for a matching xform
    found = false;
    for i=1:length(variable_transforms)
        if var1 == variable_transforms[i].from || [var1] == variable_transforms[i].from
            if var2 == variable_transforms[i].to || [var2] == variable_transforms[i].to
                found = true;
                transform_variable_values(variable_transforms[i]);
            end
        end
    end
end

function boundary(var, bid, bc_type, bc_exp=0)
    # The expression may contain variable symbols.
    # Parse it, replace variable symbols with the appropriate parts
    # Type may change in bc_exp, so convert it to an array of Any
    newbc_exp = [];
    append!(newbc_exp, bc_exp);
    if typeof(bc_exp) <: Array
        nfuns = 0;
        for i=1:length(bc_exp)
            if typeof(bc_exp[i]) == String
                ex = Meta.parse(bc_exp[i]);
                ex = replace_symbols_in_conditions(ex);
                newbc_exp[i] = string(ex);
                nfuns += makeFunctions(newbc_exp[i]);
            elseif typeof(bc_exp[i]) <: Number
                # do nothing
            else
                # What else could it be?
            end
        end
    elseif typeof(bc_exp) == String
        ex = Meta.parse(bc_exp);
        ex = replace_symbols_in_conditions(ex);
        newbc_exp = string(ex);
        nfuns = makeFunctions(newbc_exp);
    elseif typeof(bc_exp) <: Number
        nfuns = 0;
    else
        # ??
    end
    
    add_boundary_condition(var, bid, bc_type, newbc_exp, nfuns);
end

function addBoundaryID(bid, trueOnBdry)
    # trueOnBdry(x, y, z) = something # points with x,y,z on this bdry segment evaluate true here
    if typeof(trueOnBdry) == String
        trueOnBdry = stringToFunction("trueOnBdry", "x,y=0,z=0", trueOnBdry);
    end
    add_boundary_ID_to_grid(bid, trueOnBdry, grid_data);
end

function referencePoint(var, pos, val)
    add_reference_point(var, pos, val);
end

function timeInterval(T)
    prob.time_dependent = true;
    if time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    prob.end_time = T;
end

function initial(var, ics)
    nfuns = makeFunctions(ics);
    add_initial_condition(var.index, ics, nfuns);
end

function preStepFunction(fun)
    solver.set_pre_step(fun);
end
function postStepFunction(fun)
    solver.set_post_step(fun);
end

function callbackFunction(fun; name="", args=[], body="")
    if name==""
        name = string(fun);
    end
    # If args and body are not provided, this may still work internally
    # but it can't be generated for external targets.
    
    add_callback_function(CallbackFunction(name, args, body, fun));
    
    log_entry("Added callback function: "*name, 2);
end

function weakForm(var, wf)
    if typeof(var) <: Array
        # multiple simultaneous variables
        wfvars = [];
        wfex = [];
        if !(length(var) == length(wf))
            printerr("Error in weak form: # of unknowns must equal # of equations. (example: @weakform([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(wfvars, var[vi].symbol);
            push!(wfex, Meta.parse((wf)[vi]));
        end
    else
        wfex = Meta.parse(wf);
        wfvars = var.symbol;
    end
    
    log_entry("Making weak form for variable(s): "*string(wfvars));
    log_entry("Weak form, input: "*string(wf));
    
    # This is the parsing step. It goes from an Expr to arrays of Basic
    result_exprs = sp_parse(wfex, wfvars);
    if length(result_exprs) == 4 # has surface terms
        (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = result_exprs;
    else
        (lhs_symexpr, rhs_symexpr) = result_exprs;
    end
    
    # Here we set a SymExpression for each of the pieces. 
    # This is an Expr tree that is passed to the code generator.
    if length(result_exprs) == 4 # has surface terms
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, lhs_surf_symexpr, LHS, "surface");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        set_symexpressions(var, rhs_surf_symexpr, RHS, "surface");
        
        log_entry("lhs volume symexpression:\n\t"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n\t"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n\t"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n\t"*string(rhs_surf_symexpr));
        
        log_entry("Latex equation:\n\t\t \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(lhs_surf_symexpr)*
                    " ds = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(rhs_surf_symexpr)*" ds", 3);
    else
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        
        log_entry("lhs symexpression:\n\t"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n\t"*string(rhs_symexpr));
        
        log_entry("Latex equation:\n\t\t\$ \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx");
    end
    
    # change symbolic layer into code layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "volume", config.solver_type, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "volume", config.solver_type, language, gen_framework);
    if length(result_exprs) == 4
        (lhs_surf_string, lhs_surf_code) = generate_code_layer(lhs_surf_symexpr, var, LHS, "surface", config.solver_type, language, gen_framework);
        (rhs_surf_string, rhs_surf_code) = generate_code_layer(rhs_surf_symexpr, var, RHS, "surface", config.solver_type, language, gen_framework);
        log_entry("Weak form, code layer: LHS = \n"*string(lhs_string)*"\nsurfaceLHS = \n"*string(lhs_surf_string)*" \nRHS = \n"*string(rhs_string)*"\nsurfaceRHS = \n"*string(rhs_surf_string));
    else
        log_entry("Weak form, code layer: LHS = \n"*string(lhs_string)*" \n  RHS = \n"*string(rhs_string));
    end
    
    #log_entry("Julia code Expr: LHS = \n"*string(lhs_code)*" \n  RHS = \n"*string(rhs_code), 3);
    
    if language == JULIA || language == 0
        args = "args; kwargs...";
        makeFunction(args, lhs_code);
        set_lhs(var, lhs_string);
        
        makeFunction(args, rhs_code);
        set_rhs(var, rhs_string);
        
        if length(result_exprs) == 4
            makeFunction(args, string(lhs_surf_code));
            set_lhs_surface(var, lhs_surf_string);
            
            makeFunction(args, string(rhs_surf_code));
            set_rhs_surface(var, rhs_surf_string);
        end
        
    else
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

function fluxAndSource(var, fex, sex)
    flux(var, fex);
    source(var, sex);
end

function flux(var, fex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symfex = [];
        if !(length(var) == length(fex))
            printerr("Error in flux function: # of unknowns must equal # of equations. (example: flux([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symfex, Meta.parse((fex)[vi]));
        end
    else
        symfex = Meta.parse(fex);
        symvars = var.symbol;
    end
    
    log_entry("Making flux for variable(s): "*string(symvars));
    log_entry("flux, input: "*string(fex));
    
    # The parsing step
    (lhs_symexpr, rhs_symexpr) = sp_parse(symfex, symvars, is_FV=true, is_flux=true);
    
    set_symexpressions(var, lhs_symexpr, LHS, "surface");
    set_symexpressions(var, rhs_symexpr, RHS, "surface");
    
    log_entry("flux lhs symexpression:\n\t"*string(lhs_symexpr));
    log_entry("flux rhs symexpression:\n\t"*string(rhs_symexpr));
    log_entry("Latex flux equation:\n\t\t \\int_{\\partial K}"*symexpression_to_latex(lhs_symexpr)*
                    " ds - \\int_{\\partial K}"*symexpression_to_latex(rhs_symexpr)*" ds");
    
    # change symbolic layer into code layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "surface", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "surface", FV, language, gen_framework);
    log_entry("flux, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args; kwargs...";
        makeFunction(args, string(lhs_code));
        set_lhs_surface(var, lhs_string);
        
        makeFunction(args, string(rhs_code));
        set_rhs_surface(var, rhs_string);
        
    else
        set_lhs_surface(var, lhs_code);
        set_rhs_surface(var, rhs_code);
    end
end

function source(var, sex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symsex = [];
        if !(length(var) == length(sex))
            printerr("Error in source function: # of unknowns must equal # of equations. (example: source([a,b,c], [s1,s2,s3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symsex, Meta.parse((sex)[vi]));
        end
    else
        symsex = Meta.parse(sex);
        symvars = var.symbol;
    end
    
    log_entry("Making source for variable(s): "*string(symvars));
    log_entry("source, input: "*string(sex));
    
    # The parsing step
    (lhs_symexpr, rhs_symexpr) = sp_parse(symsex, symvars, is_FV=true);
    
    set_symexpressions(var, lhs_symexpr, LHS, "volume");
    set_symexpressions(var, rhs_symexpr, RHS, "volume");
    
    log_entry("source lhs symexpression:\n\t"*string(lhs_symexpr));
    log_entry("source rhs symexpression:\n\t"*string(rhs_symexpr));
    log_entry("Latex source equation:\n\t\t \\int_{K} -"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx");
    
    # change symbolic layer into code code_layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "volume", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "volume", FV, language, gen_framework);
    log_entry("source, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args; kwargs...";
        makeFunction(args, string(lhs_code));
        set_lhs(var, lhs_string);
        
        makeFunction(args, string(rhs_code));
        set_rhs(var, rhs_string);
        
    else
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

# Creates an assembly loop function for a given variable that nests the loops in the given order.
function assemblyLoops(var, indices)
    # find the associated indexer objects
    indexer_list = [];
    for i=1:length(indices)
        if typeof(indices[i]) == String
            if indices[i] == "elements" || indices[i] == "cells"
                push!(indexer_list, "elements");
            else
                for j=1:length(indexers)
                    if indices[i] == string(indexers[j].symbol)
                        push!(indexer_list, indexers[j]);
                    end
                end
            end
            
        elseif typeof(indices[i]) == Indexer
            push!(indexer_list, indices[i]);
        end
        
    end
    
    (loop_string, loop_code) = generate_assembly_loops(var, indexer_list, config.solver_type, language);
    log_entry("assembly loop function: \n\t"*string(loop_string));
    
    if language == JULIA || language == 0
        if config.solver_type == FV
            args = "var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0";
        elseif config.solver_type == CG
            args = "var, bilinear, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; rhs_only=false";
        end
        makeFunction(args, string(loop_code));
        set_assembly_loops(var, loop_string);
    else
        set_assembly_loops(var, loop_code);
    end
end

function exportCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if language == JULIA || language == 0
        file = open(filename*".jl", "w");
        println(file, "#=\nGenerated functions for "*project_name*"\n=#\n");
        for LorR in [LHS, RHS]
            if LorR == LHS
                codevol = code_strings[1];
                codesurf = code_strings[3];
            else
                codevol = code_strings[2];
                codesurf = code_strings[4];
            end
            codeassemble = code_strings[5];
            for i=1:length(variables)
                var = string(variables[i].symbol);
                if !(codevol[i] === "")
                    func_name = LorR*"_volume_function_for_"*var;
                    println(file, "function "*func_name*"(args; kwargs...)");
                    println(file, codevol[i]);
                    println(file, "end #"*func_name*"\n");
                else
                    println(file, "# No "*LorR*" volume set for "*var*"\n");
                end
                
                if !(codesurf[i] === "")
                    func_name = LorR*"_surface_function_for_"*var;
                    println(file, "function "*func_name*"(args; kwargs...)");
                    println(file, codesurf[i]);
                    println(file, "end #"*func_name*"\n");
                else
                    println(file, "# No "*LorR*" surface set for "*var*"\n");
                end
                
                if LorR == RHS
                    if !(codeassemble[i] == "" || assembly_loops[i] === nothing)
                        func_name = "assembly_function_for_"*var;
                        if config.solver_type == FV
                            args = "var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0";
                        elseif config.solver_type == CG
                            args = "var, bilinear, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; rhs_only=false";
                        end
                        println(file, "function "*func_name*"("*args*")");
                        println(file, codeassemble[i]);
                        println(file, "end #"*func_name*"\n");
                    else
                        println(file, "# No assembly function set for "*var*"\n");
                    end
                end
            end
        end
        
        close(file);
        log_entry("Exported code to "*filename*".jl");
    else
        # Should we export for other targets?
    end
end

function importCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if language == JULIA || language == 0
        file = open(filename*".jl", "r");
        lines = readlines(file, keep=true);
        for LorR in [LHS, RHS]
            if LorR == LHS
                codevol = bilinears;
                codesurf = face_bilinears;
            else
                codevol = linears;
                codesurf = face_linears;
            end
            
            # Loop over variables and check to see if a matching function is present.
            for i=1:length(variables)
                # Scan the file for a pattern like
                #   function LHS_volume_function_for_u
                #       ...
                #   end #LHS_volume_function_for_u
                #
                # Set LHS/RHS, volume/surface, and the variable name
                var = string(variables[i].symbol);
                vfunc_name = LorR*"_volume_function_for_"*var;
                vfunc_string = "";
                sfunc_name = LorR*"_surface_function_for_"*var;
                sfunc_string = "";
                afunc_name = "assembly_function_for_"*var;
                afunc_string = "";
                for st=1:length(lines)
                    if occursin("function "*vfunc_name, lines[st])
                        # s is the start of the function
                        for en=(st+1):length(lines)
                            if occursin("end #"*vfunc_name, lines[en])
                                # en is the end of the function
                                st = en; # update st
                                break;
                            else
                                vfunc_string *= lines[en];
                            end
                        end
                    elseif occursin("function "*sfunc_name, lines[st])
                        # s is the start of the function
                        for en=(st+1):length(lines)
                            if occursin("end #"*sfunc_name, lines[en])
                                # en is the end of the function
                                st = en; # update st
                                break;
                            else
                                sfunc_string *= lines[en];
                            end
                        end
                    elseif LorR == RHS && occursin("function "*afunc_name, lines[st])
                        # s is the start of the function
                        for en=(st+1):length(lines)
                            if occursin("end #"*afunc_name, lines[en])
                                # en is the end of the function
                                st = en; # update st
                                break;
                            else
                                afunc_string *= lines[en];
                            end
                        end
                    end
                end # lines loop
                
                # Generate the functions and set them in the right places
                if vfunc_string == ""
                    log_entry("Warning: While importing, no "*LorR*" volume function was found for "*var);
                else
                    makeFunction("args; kwargs...", CodeGenerator.code_string_to_expr(vfunc_string));
                    if LorR == LHS
                        set_lhs(variables[i]);
                    else
                        set_rhs(variables[i]);
                    end
                end
                if sfunc_string == ""
                    log_entry("Warning: While importing, no "*LorR*" surface function was found for "*var);
                else
                    makeFunction("args; kwargs...", CodeGenerator.code_string_to_expr(sfunc_string));
                    if LorR == LHS
                        set_lhs_surface(variables[i]);
                    else
                        set_rhs_surface(variables[i]);
                    end
                end
                if afunc_string == ""
                    log_entry("Warning: While importing, no assembly function was found for "*var*" (using default)");
                else
                    if config.solver_type == FV
                        args = "var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0";
                    elseif config.solver_type == CG
                        args = "var, bilinear, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; rhs_only=false";
                    end
                    makeFunction(args, CodeGenerator.code_string_to_expr(afunc_string));
                    set_assembly_loops(variables[i]);
                end
                
            end # vars loop
        end
        
        close(file);
        log_entry("Imported code from "*filename*".jl");
        
    else
        # TODO non-julia
    end
end

# Prints a Latex string for the equation for a variable
function printLatex(var)
    if typeof(var) <: Array
        varname = "["*string(var[1].symbol);
        for i=2:length(var)
            varname *= ", "*string(var[i].symbol);
        end
        varname *= "]";
        var = var[1];
    else
        varname = string(var.symbol);
    end
    println("Latex equation for "*varname*":");
    
    if config.solver_type == FV
        result = raw"$" * " \\frac{d}{dt}\\bar{"*string(var.symbol)*"}";
        if !(symexpressions[1][var.index] === nothing)
            result *= " + \\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
        end
        if !(symexpressions[2][var.index] === nothing)
            result *= " + \\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
        end
        result *= " = ";
        if !(symexpressions[3][var.index] === nothing)
            result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
        end
        if !(symexpressions[4][var.index] === nothing)
            if !(symexpressions[3][var.index] === nothing)
                result *= " + ";
            end
            result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
        end
        if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
            result *= "0";
        end
        result *= raw"$";
    else
        if prob.time_dependent
            result = raw"$ \left[";
            if !(symexpressions[1][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
            end
            if !(symexpressions[2][var.index] === nothing)
                if !(symexpressions[1][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
            end
            if symexpressions[1][var.index] === nothing && symexpressions[2][var.index] === nothing
                result *= "0";
            end
            result *= "\\right]_{t+dt} = \\left[";
            if !(symexpressions[3][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
            end
            if !(symexpressions[4][var.index] === nothing)
                if !(symexpressions[3][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
            end
            if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
                result *= "0";
            end
            result *= raw"\right]_{t} $";
            
        else
            result = raw"$";
            if !(symexpressions[1][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
            end
            if !(symexpressions[2][var.index] === nothing)
                if !(symexpressions[1][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
            end
            if symexpressions[1][var.index] === nothing && symexpressions[2][var.index] === nothing
                result *= "0";
            end
            result *= " = ";
            if !(symexpressions[3][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
            end
            if !(symexpressions[4][var.index] === nothing)
                if !(symexpressions[3][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
            end
            if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
                result *= "0";
            end
            result *= raw"$";
        end
    end
    println(result);
    println("");
    
    return result;
end

# Evaluate all of the initial conditions
function eval_initial_conditions()
    dim = config.dimension;

    # build initial conditions
    for vind=1:length(variables)
        if !(variables[vind].ready) && vind <= length(prob.initial)
            if !(prob.initial[vind] === nothing)
                if variables[vind].location == CELL && !(fv_grid === nothing)
                    # Need to use the fv_grid instead of grid_data
                    this_grid_data = fv_grid;
                    this_geo_factors = fv_geo_factors;
                    this_refel = fv_refel;
                else
                    this_grid_data = grid_data;
                    this_geo_factors = geo_factors;
                    this_refel = refel;
                end
                
                if typeof(prob.initial[vind]) <: Array
                    # Evaluate at nodes
                    nodal_values = zeros(length(prob.initial[vind]), size(this_grid_data.allnodes,2));
                    for ci=1:length(prob.initial[vind])
                        for ni=1:size(this_grid_data.allnodes,2)
                            if typeof(prob.initial[vind][ci]) <: Number
                                nodal_values[ci,ni] = prob.initial[vind][ci];
                            elseif dim == 1
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[ni],0,0,0);
                            elseif dim == 2
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],0,0);
                            elseif dim == 3
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],this_grid_data.allnodes[3,ni],0);
                            end
                        end
                    end
                    
                    # compute cell averages using nodal values if needed
                    if variables[vind].location == CELL
                        nel = size(this_grid_data.loc2glb, 2);
                        for ei=1:nel
                            e = elemental_order[ei];
                            glb = this_grid_data.loc2glb[:,e];
                            vol = this_geo_factors.volume[e];
                            detj = this_geo_factors.detJ[e];
                            
                            for ci=1:length(prob.initial[vind])
                                variables[vind].values[ci,e] = detj / vol * (this_refel.wg' * this_refel.Q * (nodal_values[ci,glb][:]))[1];
                            end
                        end
                    else
                        variables[vind].values = nodal_values;
                    end
                    
                else # scalar
                    nodal_values = zeros(1, size(this_grid_data.allnodes,2));
                    for ni=1:size(this_grid_data.allnodes,2)
                        if typeof(prob.initial[vind]) <: Number
                            nodal_values[ni] = prob.initial[vind];
                        elseif dim == 1
                            nodal_values[ni] = prob.initial[vind].func(this_grid_data.allnodes[ni],0,0,0);
                        elseif dim == 2
                            nodal_values[ni] = prob.initial[vind].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],0,0);
                        elseif dim == 3
                            nodal_values[ni] = prob.initial[vind].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],this_grid_data.allnodes[3,ni],0);
                        end
                    end
                    
                    # compute cell averages using nodal values if needed
                    if variables[vind].location == CELL
                        nel = size(this_grid_data.loc2glb, 2);
                        for e=1:nel
                            glb = this_grid_data.loc2glb[:,e];
                            vol = this_geo_factors.volume[e];
                            detj = this_geo_factors.detJ[e];
                            
                            variables[vind].values[e] = detj / vol * (this_refel.wg' * this_refel.Q * nodal_values[glb])[1];
                        end
                        
                    else
                        variables[vind].values = nodal_values;
                    end
                    
                end
                
                variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(variables[vind].symbol));
            end
        end
    end
end

# This will either solve the problem or generate the code for an external target.
function solve(var, nlvar=nothing; nonlinear=false)
    if use_cachesim
        return cachesimSolve(var);
    end
    
    global time_stepper; # This should not be necessary. It will go away eventually
    
    # Generate files or solve directly
    if !(!generate_external && (language == JULIA || language == 0)) # if an external code gen target is ready
        if typeof(var) <: Array
            varind = var[1].index;
        else
            varind = var.index;
        end
        generate_all_files(var, bilinears[varind], face_bilinears[varind], linears[varind], face_linears[varind]);
        
    else
        if typeof(var) <: Array
            varnames = "["*string(var[1].symbol);
            for vi=2:length(var)
                varnames = varnames*", "*string(var[vi].symbol);
            end
            varnames = varnames*"]";
            varind = var[1].index;
        else
            varnames = string(var.symbol);
            varind = var.index;
        end
        
        # Evaluate initial conditions if not already done
        eval_initial_conditions();
        
        # If any of the variables are indexed types, generate assembly loops
        var_indexers = [];
        if typeof(var) <: Array
            for vi=1:length(var)
                if !(var[vi].indexer === nothing)
                    push!(var_indexers, var[vi].indexer);
                end
            end
        else
            if !(var.indexer === nothing)
                var_indexers = [var.indexer];
            end
        end
        if length(var_indexers)>0 && assembly_loops[varind] === nothing
            log_entry("Indexed variables detected, but no assembly loops specified. Using default.")
            assemblyLoops(var, ["elements"; var_indexers]);
        end
        
        # Use the appropriate solver
        if config.solver_type == CG
            lhs = bilinears[varind];
            rhs = linears[varind];
            loop_func = assembly_loops[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
                if (nonlinear)
                    if time_stepper.type == EULER_EXPLICIT || time_stepper.type == LSRK4
                        printerr("Warning: Use implicit stepper for nonlinear problem. cancelling solve.");
                        return;
                    end
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs, time_stepper, assemble_func=loop_func));
				else
                	t = @elapsed(result = CGSolver.linear_solve(var, lhs, rhs, time_stepper, assemble_func=loop_func));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs, assemble_func=loop_func));
                else
                    t = @elapsed(result = CGSolver.linear_solve(var, lhs, rhs, assemble_func=loop_func));
				end
                
                # place the values in the variable value arrays
                if typeof(var) <: Array && length(result) > 1
                    tmp = 0;
                    totalcomponents = 0;
                    for vi=1:length(var)
                        totalcomponents = totalcomponents + length(var[vi].symvar);
                    end
                    for vi=1:length(var)
                        components = length(var[vi].symvar);
                        for compi=1:components
                            #println("putting result "*string(compi+tmp)*":"*string(totalcomponents)*":end in var["*string(vi)*"].values["*string(compi)*"]")
                            var[vi].values[compi,:] = result[(compi+tmp):totalcomponents:end];
                        end
                        tmp = tmp + components;
                    end
                elseif length(result) > 1
                    components = length(var.symvar);
                    for compi=1:components
                        var.values[compi,:] = result[compi:components:end];
                    end
                end
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == DG
            lhs = bilinears[varind];
            rhs = linears[varind];
            slhs = face_bilinears[varind];
            srhs = face_linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	# t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs, time_stepper));
                    printerr("Nonlinear solver not ready for DG");
                    return;
				else
                	t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs, time_stepper));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs));
                else
                    t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs));
				end
                
                # place the values in the variable value arrays
                if typeof(var) <: Array && length(result) > 1
                    tmp = 0;
                    totalcomponents = 0;
                    for vi=1:length(var)
                        totalcomponents = totalcomponents + length(var[vi].symvar);
                    end
                    for vi=1:length(var)
                        components = length(var[vi].symvar);
                        for compi=1:components
                            var[vi].values[compi,:] = result[(compi+tmp):totalcomponents:end];
                            tmp = tmp + 1;
                        end
                    end
                elseif length(result) > 1
                    components = length(var.symvar);
                    for compi=1:components
                        var.values[compi,:] = result[compi:components:end];
                    end
                end
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == FV
            slhs = bilinears[varind];
            srhs = linears[varind];
            flhs = face_bilinears[varind];
            frhs = face_linears[varind];
            loop_func = assembly_loops[varind];
            
            if prob.time_dependent
                if fv_grid === nothing
                    time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                else
                    time_stepper = init_stepper(fv_grid.allnodes, time_stepper);
                end
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	printerr("Nonlinear solver not ready for FV");
                    return;
				else
                	t = @elapsed(result = solver.linear_solve(var, slhs, srhs, flhs, frhs, time_stepper, loop_func));
				end
                # result is already stored in variables
                log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            else
                # does this make sense?
                printerr("FV assumes time dependence. Set initial conditions etc.");
            end
        end
    end
    
    # At this point all of the final values for the variables in var should be placed in var.values.
end

# When using cachesim, this will be used to simulate the solve.
function cachesimSolve(var, nlvar=nothing; nonlinear=false)
    if !(!generate_external && (language == JULIA || language == 0))
        printerr("Cachesim solve is only ready for Julia direct solve");
    else
        if typeof(var) <: Array
            varnames = "["*string(var[1].symbol);
            for vi=2:length(var)
                varnames = varnames*", "*string(var[vi].symbol);
            end
            varnames = varnames*"]";
            varind = var[1].index;
        else
            varnames = string(var.symbol);
            varind = var.index;
        end

        lhs = bilinears[varind];
        rhs = linears[varind];
        
        t = @elapsed(result = linear_solve_cachesim(var, lhs, rhs));
        log_entry("Generated cachesim ouput for "*varnames*".(took "*string(t)*" seconds)", 1);
    end
end

# Writes the values for each variable to filename in the given format.
function output_values(vars, filename; format="raw", ascii=false)
    available_formats = ["raw", "csv", "vtk"];
    if format == "vtk"
        output_values_vtk(vars, filename, ascii)
        log_entry("Writing values to file: "*filename*".vtu");
        
    else
        filename *= "."*format;
        file = open(filename, "w");
        log_entry("Writing values to file: "*filename);
        
        if format == "raw"
            output_values_raw(vars, file)
        elseif format == "csv"
            output_values_csv(vars, file)
        else
            println("Unknown output file format("*string(format)*"). Choose from: "*string(available_formats));
        end
        log_entry("Finished writing to "*filename);
        
        close(file);
    end
end

function finalize_finch()
    # Finalize generation
    finalize_code_generator();
    if use_cachesim
        CachesimOut.finalize();
    end
    # log
    close_log();
    # # mpi
    # if config.use_mpi
    #     MPI.Finalize(); # If MPI is finalized, it cannot be reinitialized without restarting Julia
    # end
    if config.proc_rank == 0
        println("Finch has completed.");
    end
end

### Other specialized functions ###

function cachesim(use)
    log_entry("Using cachesim - Only cachesim output will be generated.", 1);
    global use_cachesim = use;
end

function morton_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive!(grid_data, griddim, MORTON_ORDERING));
    log_entry("Reordered nodes to Morton. Took "*string(t)*" sec.", 2);
end

function morton_elements(griddim)
    global elemental_order = get_recursive_order(MORTON_ORDERING, config.dimension, griddim);
    log_entry("Reordered elements to Morton.", 2);
    ef_nodes();
end

function hilbert_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive!(grid_data, griddim, HILBERT_ORDERING));
    log_entry("Reordered nodes to Hilbert. Took "*string(t)*" sec.", 2);
end

function hilbert_elements(griddim)
    global elemental_order = get_recursive_order(HILBERT_ORDERING, config.dimension, griddim);
    log_entry("Reordered elements to Hilbert.", 2);
    ef_nodes();
end

function tiled_nodes(griddim, tiledim)
    t = @elapsed(global grid_data = reorder_grid_tiled(grid_data, griddim, tiledim));
    log_entry("Reordered nodes to tiled. Took "*string(t)*" sec.", 2);
end

function tiled_elements(griddim, tiledim)
    global elemental_order = get_tiled_order(config.dimension, griddim, tiledim, true);
    log_entry("Reordered elements to tiled("*string(tiledim)*").", 2);
    ef_nodes();
end

function ef_nodes()
    t = @elapsed(global grid_data = reorder_grid_element_first!(grid_data, config.basis_order_min, elemental_order));
    log_entry("Reordered nodes to EF. Took "*string(t)*" sec.", 2);
end

function random_nodes(seed = 17)
    t = @elapsed(global grid_data = reorder_grid_random!(grid_data, seed));
    log_entry("Reordered nodes to random. Took "*string(t)*" sec.", 2);
end

function random_elements(seed = 17)
    global elemental_order = random_order(size(grid_data.loc2glb,2), seed);
    log_entry("Reordered elements to random.", 2);
    random_nodes(seed);
end