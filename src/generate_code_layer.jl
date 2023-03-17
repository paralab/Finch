#=
Use the IR to generate code
Redirects to solver and target specific functions.
=#

function generate_code_layer(var, IR, solver, language, framework)
    # if use_cachesim
    #     return generate_code_layer_cachesim(var, entities, terms, lorr, vors);
    # end
    
    if language == JULIA
        if finch_state.config.use_gpu
            (code, aux_code) = generate_code_layer_julia_gpu(var, IR, solver);
            
        else
            # code holds the solve function
            # aux_code holds any other code that needs to be imported for this to work
            (code, aux_code) = generate_code_layer_julia(var, IR, solver);
        end
        
    ### External targets ##############################################################
    # The appropriate code gen function should be set
    else
        code = code_gen_context.external_generate_code_files_function(var, IR);
        # We won't actually return code because it has been written to files.
        code = fill("# See code in generated files ", length(var));
        aux_code = "";
    end
    
    return (code, aux_code);
end
