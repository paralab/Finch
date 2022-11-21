#=
Use the IR to generate code
Redirects to solver and target specific functions.
=#

function generate_code_layer(var, IR, solver, language, framework)
    # if use_cachesim
    #     return generate_code_layer_cachesim(var, entities, terms, lorr, vors);
    # end
    
    if language == JULIA
        code = generate_code_layer_julia(var, IR, solver);
        
    ### External targets ##############################################################
    # The appropriate code gen function should be set
    else
        code = code_gen_context.external_generate_code_files_function(var, IR);
        # We won't actually return code because it has been written to files.
        code = "# See code in generated files ";
    end
    
    return code;
end
