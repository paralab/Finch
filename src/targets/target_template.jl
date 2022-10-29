#=
This contains all of the pieces needed to add a new code gen target,
which generates a set of code files from a given IR.

The following two functions are required for a target:
1. get_external_language_elements() - file extensions, comment chars etc.
2. generate_external_files(var, IR) - Writes all files based on the configuration and IR

You will have access to anything in the Finch module scope because this file will be included there.
=#

# Returns a set of basic language elements such as comment characters and file extensions.
function get_external_language_elements()
    file_extension = ".jl";
    comment_char = "#";
    block_comment = ["#="; "=#"];
    
    return (file_extension, comment_char, block_comment);
end

# Generates all of the code files from the given IR.
# This will most likely split up the task among many other functions.
function generate_external_files(var, IR) 
    
end