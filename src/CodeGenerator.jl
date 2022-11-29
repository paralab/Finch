#=
Module for code generation
=#
module CodeGenerator

export init_code_generator, finalize_code_generator, set_generation_target,
        generate_all_files, add_generated_file,
        # generate_main, generate_config, generate_prob, generate_mesh, generate_genfunction, 
        # generate_bilinear, generate_linear, generate_stepper, generate_output,
        generate_code_layer, generate_assembly_loops
        #, generate_code_layer_surface, generate_code_layer_fv

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

# IR symbols
import ..Finch: IR_entry_types, IR_string, print_IR, repr_IR,
            IR_part, IR_data_node, IR_data_access, IR_operation_node, 
            IR_block_node, IR_loop_node, IR_conditional_node, IR_comment_node
import ..Finch: generate_linalg_TDM_product, generate_the_one_pattern, generate_linalg_Tv_product,
            apply_indexed_access, generate_difference_norms_and_update, generate_residual_norms_and_update, 
            generate_local_to_global_fem
#
# Temporary placeholders for external code gen functions that must be provided.
# These are reassigned in set_custom_target()
function default_language_elements_function() return (".jl", "#", ["#=", "=#"]) end;
function default_code_files_function() end;
                
# Holds basic generator info
mutable struct CodeGenContext
    gen_dir::String
    gen_file_name::String
    file_extension::String
    comment_char::String
    block_comment_char::Vector{String}
    header_text::String
    
    gen_files::Vector
    external_get_language_elements_function::Function
    external_generate_code_files_function::Function
    
    using_custom_target::Bool
    
    CodeGenContext() = new(
        "", "", "", "", ["",""], "", [], 
        default_language_elements_function, default_code_files_function, false
    )
end

# general code generator functions
include("code_generator_utils.jl");
include("generate_code_layer.jl");

# code gen functions for each solver type and target
include("generate_code_layer_julia.jl");

# # target specific code gen functions
# include("generate_code_layer_dendro.jl");
# include("generate_code_layer_homg.jl");
# include("generate_code_layer_matlab.jl");
# include("generate_code_layer_cachesim.jl");

#Matlab
# include("generate_matlab_utils.jl");
# include("generate_matlab_files.jl");
# include("generate_homg_files.jl");
# #C++
# include("generate_cpp_utils.jl");
# include("generate_dendro_files.jl");


#### Note
# default Dendro parameters
# parameters = (5, 1, 0.3, 0.000001, 100);#(maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
####

# a global context
code_gen_context = CodeGenContext();

function init_code_generator(dir, name, header)
    code_gen_context.gen_dir = dir;
    code_gen_context.gen_file_name = name;
    code_gen_context.header_text = header;
end

# Sets the functions to be used during external code generation
function set_generation_target(lang_elements, file_maker)
    code_gen_context.external_get_language_elements_function = lang_elements;
    code_gen_context.external_generate_code_files_function = file_maker;
    code_gen_context.using_custom_target = true;
    
    (file_extension, comment_char, block_comment_char) = Base.invokelatest(code_gen_context.external_get_language_elements_function);
    code_gen_context.file_extension = file_extension;
    code_gen_context.comment_char = comment_char;
    code_gen_context.block_comment_char = block_comment_char;
end

function add_generated_file(filename; dir="", make_header_text=true)
    if length(dir) > 0
        code_dir = code_gen_context.gen_dir*"/"*dir;
        if !isdir(code_dir)
            mkdir(code_dir);
        end
    else
        code_dir = code_gen_context.gen_dir;
    end
    newfile = open(code_dir*"/"*filename, "w");
    push!(code_gen_context.gen_files, newfile);
    if make_header_text
        generate_head(newfile, code_gen_context.header_text);
    end
    
    return newfile;
end

function generate_all_files(var, IR; parameters=0)
    if code_gen_context.using_custom_target
        code_gen_context.external_generate_code_files_function(var, IR);
    end
end

function finalize_code_generator()
    for f in code_gen_context.gen_files
        close(f);
    end
    log_entry("Closed generated code files.");
end

#### Utilities ####

function comment(file,line)
    println(file, code_gen_context.comment_char * line);
end

function commentBlock(file,text)
    print(file, "\n"*code_gen_context.block_comment_char[1]*"\n"*text*"\n"*code_gen_context.block_comment_char[2]*"\n");
end

function generate_head(file, text)
    comment(file,"This file was generated by Finch.");
    commentBlock(file, text);
end

# for writing structs to binary files
# format is | number of structs[Int64] | sizes of structs[Int64*num] | structs |
function write_binary_head(f, num, szs)
    Nbytes = 0;
    write(f, num);
    Nbytes += sizeof(num)
    for i=1:length(szs)
        write(f, szs[i])
        Nbytes += sizeof(szs[i])
    end
    return Nbytes;
end

# Write an array to a binary file.
# Return number of bytes written.
function write_binary_array(f, a, with_counts=false, zero_index=false)
    Nbytes = 0;
    if with_counts
        write(f, Int64(length(a)));
        if length(a) > 0 && isbits(a[1])
            write(f, Int64(sizeof(a[1])));
        else # empty aray or array of arrays has element size = 0
            write(f, Int64(0));
        end
    end
    for i=1:length(a)
        if isbits(a[i])
            if zero_index && typeof(a[i]) <: Integer # This is probably an index. Change from 1-based to 0-based.
                write(f, a[i]-1);
            else
                write(f, a[i]);
            end
            Nbytes += sizeof(a[i]);
        else
            Nbytes += write_binary_array(f,a[i], with_counts, zero_index);
        end
    end
    return Nbytes;
end

# Assumes that the struct only has isbits->true types or arrays.
# Returns number of bytes written.
# with_counts=true will add number of pieces and size of pieces before each piece.(Int64, Int64)
# zero_index=true will attempt to change indices from 1-based to 0-based
function write_binary_struct(f, s, with_counts=false, zero_index=false)
    Nbytes = 0;
    for fn in fieldnames(typeof(s))
        comp = getfield(s, fn);
        if isbits(comp)
            if with_counts
                write(f, Int64(1));
                write(f, Int64(sizeof(comp)));
            end
            write(f, comp);
            Nbytes += sizeof(comp)
        else
            Nbytes += write_binary_array(f,comp, with_counts, zero_index);
        end
    end
    return Nbytes;
end

# Write the grid to a binary file intended to be imported in C++
# This includes various extra numbers to size the arrays.
function write_grid_to_file(file, grid)
    # various numbers
    write(file, Int(size(grid.allnodes,1)));    # dimension
    write(file, Int(size(grid.loc2glb,2)));     # This is local nel (owned + ghost)
    write(file, Int(size(grid.allnodes,2)));    # nnodes local
    
    write(file, grid.nel_global);               # These are written even for non-partitioned grid
    write(file, grid.nnodes_global);            #
    
    write(file, Int(size(grid.loc2glb,1)));     # nodes per element
    write(file, Int(size(grid.glbvertex,1)));   # vertices per element
    write(file, Int(size(grid.element2face,1)));# faces per element
    
    write(file, Int(size(grid.face2element,2)));# nfaces
    write(file, Int(size(grid.face2glb,1)));    # nodes per face
    
    # Now the data in grid
    write_binary_array(file, grid.allnodes, true);
    write_binary_array(file, grid.bdry, true, true);
    write_binary_array(file, grid.bdryface, true, true);
    write_binary_array(file, grid.bdrynorm, true);
    write_binary_array(file, grid.bids, true);
    write_binary_array(file, grid.loc2glb, true, true);
    write_binary_array(file, grid.glbvertex, true, true);
    write_binary_array(file, grid.face2glb, true, true);
    write_binary_array(file, grid.element2face, true, true);
    write_binary_array(file, grid.face2element, true, true);
    write_binary_array(file, grid.facenormals, true);
    write_binary_array(file, grid.faceRefelInd, true, true);
    write_binary_array(file, grid.facebid, true);
    
    write(file, Int(grid.is_subgrid));
    write(file, grid.nel_owned);
    write(file, grid.nel_ghost);
    write(file, grid.nface_owned);
    write(file, grid.nface_ghost);
    
    write(file, grid.nnodes_borrowed);
    write_binary_array(file, grid.grid2mesh, true);
    
    if grid.nel_ghost == 0 # FE only
        write_binary_array(file, grid.partition2global, true, true);
        write_binary_array(file, grid.node_owner, true);
        write_binary_array(file, grid.global_bdry_index, true);
    end
    
    if grid.nel_ghost > 0 # FV only
        write_binary_array(file, grid.element_owner, true);
        write(file, grid.num_neighbor_partitions);
        write_binary_array(file, grid.neighboring_partitions, true);
        write_binary_array(file, grid.ghost_counts, true);
        write_binary_array(file, grid.ghost_index, true);
    end
    
end

# Write the refel to a binary file intended to be imported in C++
function write_refel_to_file(file, refel)
    write(file, refel.dim);     # dimension
    write(file, refel.N);       # order
    write(file, refel.Np);      # number of nodes
    write(file, refel.Nqp);     # number of quadrature points
    write(file, refel.Nfaces);  # number of faces
    write_binary_array(file, refel.Nfp, true); # number of face points per face
    # nodes and vandermonde
    write_binary_array(file, refel.r, true);
    write_binary_array(file, refel.wr, true);
    write_binary_array(file, refel.g, true);
    write_binary_array(file, refel.wg, true);
    write_binary_array(file, Array(transpose(refel.V)), true); ### These need to be transposed  for row-major -> col-major
    write_binary_array(file, Array(transpose(refel.gradV)), true);
    write_binary_array(file, Array(transpose(refel.invV)), true);
    write_binary_array(file, Array(transpose(refel.Vg)), true);
    write_binary_array(file, Array(transpose(refel.gradVg)), true);
    write_binary_array(file, Array(transpose(refel.invVg)), true);
    # quadrature matrices
    write_binary_array(file, Array(transpose(refel.Q)), true);
    write_binary_array(file, Array(transpose(refel.Qr)), true);
    write_binary_array(file, Array(transpose(refel.Qs)), true);
    write_binary_array(file, Array(transpose(refel.Qt)), true);
    write_binary_array(file, Array(transpose(refel.Ddr)), true);
    write_binary_array(file, Array(transpose(refel.Dds)), true);
    write_binary_array(file, Array(transpose(refel.Ddt)), true);
    
    # surface versions
    write_binary_array(file, refel.face2local, true, true);
    write_binary_array(file, refel.surf_r, true);
    write_binary_array(file, refel.surf_wr, true);
    write_binary_array(file, refel.surf_g, true);
    write_binary_array(file, refel.surf_wg, true);
    write_binary_array(file, Array(transpose(refel.surf_V)), true);
    write_binary_array(file, Array(transpose(refel.surf_gradV)), true);
    write_binary_array(file, Array(transpose(refel.surf_Vg)), true);
    write_binary_array(file, Array(transpose(refel.surf_gradVg)), true);
    
    write_binary_array(file, Array(transpose(refel.surf_Q)), true);
    write_binary_array(file, Array(transpose(refel.surf_Qr)), true);
    write_binary_array(file, Array(transpose(refel.surf_Qs)), true);
    write_binary_array(file, Array(transpose(refel.surf_Qt)), true);
    write_binary_array(file, Array(transpose(refel.surf_Ddr)), true);
    write_binary_array(file, Array(transpose(refel.surf_Dds)), true);
    write_binary_array(file, Array(transpose(refel.surf_Ddt)), true);
end

# Write the geometric factors to a binary file intended to be imported in C++
function write_geometric_factors_to_file(file, geofacs)
    if size(geofacs.detJ,1) > 1
        write(file, Int8(1)); # Constant jacobian
    else
        write(file, Int8(0)); # NOT Constant jacobian
    end
    
    # number of elements
    write(file, Int(size(geofacs.detJ, 2)));
    # number of values per element
    write(file, Int(size(geofacs.detJ, 1)));
    
    write_binary_array(file, geofacs.detJ, true);
    
    for i=1:length(geofacs.J)
        write_binary_array(file, geofacs.J[i].rx, true);
        write_binary_array(file, geofacs.J[i].ry, true);
        write_binary_array(file, geofacs.J[i].rz, true);
        write_binary_array(file, geofacs.J[i].sx, true);
        write_binary_array(file, geofacs.J[i].sy, true);
        write_binary_array(file, geofacs.J[i].sz, true);
        write_binary_array(file, geofacs.J[i].tx, true);
        write_binary_array(file, geofacs.J[i].ty, true);
        write_binary_array(file, geofacs.J[i].tz, true);
    end
    
    write_binary_array(file, geofacs.volume, true);
    write_binary_array(file, geofacs.area, true);
    write_binary_array(file, geofacs.face_detJ, true);
end

end # module