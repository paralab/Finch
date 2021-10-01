#=
This contains all of the pieces needed to add a new code gen target.
The following three functions must be provided.

1. get_external_language_elements() - file extensions, comment chars etc.
2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code

You will have access to anything in the Finch module scope because this file will be included there.
=#

function get_external_language_elements()
    file_extension = ".jl";
    comment_char = "#";
    block_comment = ["#="; "=#"];
    # I imagine there could be more here in the future. Maybe not.
    
    return (file_extension, comment_char, block_comment);
end

#=
Translate the symbolic layer into code for the elemental assembly.
Input:
    var - unknown variable or array of variables.
    entities - leaf nodes of the computational graph. Type is SymEntity(see symexpression.jl) or number.
    terms - array of additive terms, separated for convenience. [t1, t2, t3] for expression t1+t2+t3 (see note below)
    lorr - LHS or RHS will have constant values LHS or RHS which are "lhs", "rhs"
    vors - volume or surface will have values "volume" or "surface"

Note on terms: Each term is a SymExpression(see symexpression.jl) which is essentially a Julia Expr and list of entities.
The Expr holds the computational graph for each term. LHS and RHS are determined by separating terms that involve 
unknowns from those that don't. Also, in well formed FEM problems each term should be multiplied by a test function,
and thus quadrature can be done for each term separately.
=#
function generate_external_code_layer(var, entities, terms, lorr, vors)
    # The whole body of code is stored in a string that will be inserted in the 
    # appropriate function during file writing.
    code = "";
    
    # ... code generation here ...
    
    return code;
end

# Create all of the code files.
# The input are code strings generated above. They will be inserted into their appropriate files.
function generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf)
    # This can be split up into smaller functions as you wish.
    # generate_custom_utils_file();
    # generate_custom_main_file();
    # generate_custom_config_file();
    # etc.
end
