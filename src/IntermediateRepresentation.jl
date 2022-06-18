#=
The intermediate representation.
It should:
- Give code generators enough info to generate code without depending
  on DSL level info.
- Be independent of code target.

It is a tree (implements AbstractTrees)

There are different types of nodes:
- data - some piece of data with a name and optional index
- operation - allocation, evaluation, assignment, etc
- block - a list of IR_parts grouped together
- loop - a loop containing a block
- conditional - a conditional containing a block and optional else block
=#

module IntermediateRepresentation

# Implement AbstractTrees for plotting/printing/using other tools if desired
using AbstractTrees

export IR_entry_types, IR_string, print_tree, testIR
export IR_part, IR_data_node, IR_data_access, IR_operation_node, IR_block_node, IR_loop_node, IR_conditional_node
export build_IR_fem

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

AbstractTrees.children(a::Nothing) = ();
AbstractTrees.children(a::Symbol) = ();
AbstractTrees.children(a::Number) = ();
AbstractTrees.printnode(io::IO, a::Nothing) = print(io,"nothing");
AbstractTrees.printnode(io::IO, a::Symbol) = print(io,string(a));
AbstractTrees.printnode(io::IO, a::Number) = print(io,string(a));

# All of the node types are subtypes of this.
abstract type IR_part end;

# This acts as an enum for IR entry types. It is constructed or passed where needed
# rather than using @enum, which creates a bunch of global names.
# Create with `types = IR_entry_types();`
struct IR_entry_types
    name::Dict{Int8, String} # A way to get the name of an entry type from its value
    
    # loop types
    time_loop::Int8   # = 1 # time stepping
    space_loop::Int8  # = 2 # element/face/node
    dof_loop::Int8    # = 3 # variable components
    index_loop::Int8  # = 4 # other indices such as matrix indices
    
    # operation types
    allocate_op::Int8 # = 11 # allocation has at least two args: type, dims...
    assign_op::Int8     # = 12 # assignment is like lhs = rhs, lhs is a data node, rhs is data or op
    function_op::Int8   # = 13 # function evaluation has name and tuple of args like name(args)
    math_op::Int8       # = 14 # a function_op for arithmetic just in case it is useful to specify
    named_op::Int8      # = 15 # a special op keyword that will be interpreted by codegen as needed
    
    # data types
    scalar_data::Int8    # = 21 # a single value (x)
    array_data::Int8     # = 22 # an array location (x[1,3])
    matrix_data::Int8    # = 23 # an entire matrix (A) size is stored in index
    vector_data::Int8    # = 24 # an entire vector (v) size is stored in index
    int_32_data::Int8    # = 25
    int_64_data::Int8    # = 26
    float_32_data::Int8  # = 27
    float_64_data::Int8  # = 28
    
    # access types
    read_access::Int8  # = 31 
    write_access::Int8 # = 32
    
    IR_entry_types() = new(
        Dict{Int8, String}(
            1=>"time loop",
            2=>"space loop",
            3=>"dof loop",
            4=>"index loop",
            
            11=>"allocate op",
            12=>"assign op",
            13=>"function op",
            14=>"math op",
            15=>"named op",
            
            21=>"scalar data",
            22=>"array data",
            23=>"matrix data",
            24=>"vector_data",
            25=>"Int32",
            26=>"Int64",
            27=>"Float32",
            28=>"Float64",
            
            31=>"read",
            32=>"write"
        ),
        1,2,3,4, 11,12,13,14,15, 21,22,23,24,25,26,27,28, 31,32
        );
end

# This represents a location for data.
# It could be a variable name or an array name with index symbols.
struct IR_data_node <: IR_part
    type::Int8 # The type of data
    var::Symbol
    index::Vector{Union{Symbol,Int,IR_part}}
    IR_data_node(t::Int8, v::Symbol) = new(t,v,[]);
    IR_data_node(t::Int8, v::Symbol, i::Vector) = new(t,v,i);
end
AbstractTrees.children(a::IR_data_node) = a.index;
AbstractTrees.printnode(io::IO, a::IR_data_node) = print(io, string(a.var));

# Holds info about a data access
# type is read or write
# data is an IR_data_node which has variable symbol and index
struct IR_data_access <: IR_part
    type::Int8 # The type of access
    data::IR_data_node
end
AbstractTrees.children(a::IR_data_access) = ();
AbstractTrees.printnode(io::IO, a::IR_data_access) = AbstractTrees.printnode(io, a.data);

# An operation could be an allocation, assignment, function eval, or special named op
# It has a type and a vector of args. The arg pattern will depend on the type usually.
struct IR_operation_node <: IR_part
    type::Int8   # The type of operation
    args::Vector #{Union{IR_part,Symbol,Number,Nothing}}
    deps::Vector{IR_data_access}
    
    IR_operation_node(t::Int8, a::Vector) =  #{Union{IR_part,Symbol,Number,Nothing}}) = 
        new(t, a, get_accesses(a));
end
function AbstractTrees.children(a::IR_operation_node)
    IRtypes = IR_entry_types();
    if a.type == IRtypes.allocate_op
        return a.args[2:end];
    elseif a.type == IRtypes.assign_op
        return a.args;
    elseif a.type == IRtypes.function_op
        return a.args[2:end];
    elseif a.type == IRtypes.named_op
        return a.args[2:end];
    else
        return a.args
    end
end
function AbstractTrees.printnode(io::IO, a::IR_operation_node)
    IRtypes = IR_entry_types();
    if a.type == IRtypes.allocate_op
        print(io, "Allocate{"*IRtypes.name[a.args[1]]*"}");
    elseif a.type == IRtypes.assign_op
        print(io,"=");
    elseif a.type == IRtypes.function_op
        print(io, IR_string(a.args[1]));
    elseif a.type == IRtypes.named_op
        print(io, IR_string(a.args[1]));
    else
        print(io, "unknown op");
    end
end

# This is a block of IR_parts such as a set of statements.
# Use for the body of a loop or conditional or for a top-level container
mutable struct IR_block_node <: IR_part
    parts::Vector{IR_part}
end
AbstractTrees.children(a::IR_block_node) = a.parts;
AbstractTrees.printnode(io::IO, a::IR_block_node) = print(io,"block");

# This is a loop.
# The type could be time, space, dof, or index
# It has a name of the collection it is looping over, the iterator symbol, and limits.
# If limits are both zero it loops over the whole collection.
mutable struct IR_loop_node <: IR_part
    type::Int8   # The type of loop
    collection::Symbol
    iterator::Symbol
    first::Int
    last::Int
    body::IR_block_node #Vector{IR_part}
    deps::Vector{IR_data_access}
    function IR_loop_node(t::Int8, c::Symbol, i::Symbol, f::Int, l::Int, b)
        if typeof(b) == IR_block_node
            block = b;
        elseif typeof(b) <: Vector
            block = IR_block_node(b);
        else
            block = IR_block_node([b]);
        end
        
        return new(t,c,i,f,l,block,get_accesses(b));
    end
end
AbstractTrees.children(a::IR_loop_node) = (a.body,);
function AbstractTrees.printnode(io::IO, a::IR_loop_node)
    IRtypes = IR_entry_types();
    str = "Loop{";
    if a.first == 0 && a.last == 0
        str *= string(a.collection);
    else
        str *= string(a.iterator)*"="*string(a.first)*", "*string(a.last);
    end
    str *= "}";
    print(io, str);
end

# This is a conditional block
# It has a condition, true-body, else-body(optional)
mutable struct IR_conditional_node <: IR_part
    condition::Union{IR_data_node, IR_operation_node}
    body::IR_block_node #Vector{IR_part}
    elsepart::Union{IR_block_node, Nothing}
    deps::Vector{IR_data_access}
    function IR_conditional_node(c, b, e)
        if typeof(b) == IR_block_node
            block = b;
        elseif typeof(b) <: Vector
            block = IR_block_node(b);
        else
            block = IR_block_node([b]);
        end
        if typeof(e) == IR_block_node || typeof(e) == Nothing
            elseblock = e;
        elseif typeof(e) <: Vector
            elseblock = IR_block_node(e);
        else
            elseblock = IR_block_node([e]);
        end
        new(c,block,elseblock,append!(get_accesses(b), get_accesses(e)));
    end
end
AbstractTrees.children(a::IR_conditional_node) = (a.condition, a.body, a.elsepart);
AbstractTrees.printnode(io::IO, a::IR_conditional_node) = a.elsepart===nothing ? print(io, "if") : print(io, "if/else");

function get_accesses(part, is_write=false)
    access::Vector{IR_data_access} = [];
    IRtypes = IR_entry_types();
    
    if typeof(part) <: IR_block_node
        for p in part.parts
            append!(access, get_accesses(p));
        end
        
    elseif typeof(part) <: Array
        for p in part
            append!(access, get_accesses(p));
        end
        
    elseif typeof(part) == IR_data_node
        if is_write
            access = [IR_data_access(IRtypes.write_access, part)];
        else
            access = [IR_data_access(IRtypes.read_access, part)];
        end
        
    elseif typeof(part) == IR_conditional_node || typeof(part) == IR_loop_node
        append!(access, part.deps);
        
    elseif typeof(part) == IR_operation_node
        append!(access, part.deps);
    end
    
    return access;
end

Base.show(io::IO, x::IR_part) = print(io, IR_string(x));
Base.string(x::IR_part) = IR_string(x);
function IR_string(a::IR_block_node)
    result = "";
    for b in a.parts
        bodylines = split(string(b), "\n");
        for s in bodylines
            result *= "    "*s*"\n";
        end
    end
    return result;
end
function IR_string(a::IR_loop_node)
    if a.first==0 && a.last==0
        result = "for "*string(a.iterator)*" = "*string(a.collection)*"\n";
    else
        result = "for "*string(a.iterator)*" = "*string(a.first)*":"*string(a.last)*"\n";
    end
    # for b in a.body
    #     bodylines = split(string(b), "\n");
    #     for s in bodylines
    #         result *= "    "*s*"\n";
    #     end
    # end
    result *= string(a.body);
    result *= "end\n";
    return result;
end
function IR_string(a::IR_conditional_node)
    result = "if "*string(a.condition)*"\n";
    # for i=1:length(a.body)
    #     bodylines = split(string(a.body[i]), "\n");
    #     for s in bodylines
    #         result *= "    "*s*"\n";
    #     end
    # end
    result *= string(a.body);
    if !(a.elsepart === nothing) && length(a.parts) > 0
        result *= "else";
        if typeof(a.elsepart.parts[1]) != IR_conditional_node
            result *= "\n";
        end
        # for i=1:length(a.elsepart)
        #     bodylines = split(string(a.elsepart[i]), "\n");
        #     for s in bodylines
        #         result *= "    "*s*"\n";
        #     end
        # end
        result *= string(a.elsepart);
    end
    result *= "end\n";
    return result;
end
function IR_string(a::IR_operation_node)
    IRtypes = IR_entry_types();
    if a.type == IRtypes.allocate_op
        if a.args[1] == IRtypes.int_32_data
            typename = "Int32";
        elseif a.args[1] == IRtypes.int_64_data
            typename = "Int64";
        elseif a.args[1] == IRtypes.float_32_data
            typename = "Float32";
        elseif a.args[1] == IRtypes.float_64_data
            typename = "Float64";
        else
            typename = "UNKNOWNTYPE";
        end
        dimstring = string(a.args[2]);
        for i=3:length(a.args)
            dimstring *= ", "*string(a.args[i]);
        end
        result = "zeros("*typename*", "*dimstring*")";
        
    elseif a.type == IRtypes.assign_op
        result = string(a.args[1])*" = "*string(a.args[2])*";\n";
        
    elseif a.type == IRtypes.function_op || a.type == IRtypes.math_op || a.type == IRtypes.named_op
        result = string(a.args[1])*"(";
        for i = 2:length(a.args)
            result *= string(a.args[i])
            if i < length(a.args) result *= ", "; end
        end
        result *= ")";
        
    else
        result = "";
    end
    
    return result;
end
function IR_string(a::IR_data_node)
    IRtypes = IR_entry_types();
    result = string(a.var);
    if a.type==IRtypes.array_data && length(a.index) > 0
        result *= "["*string(a.index[1]);
        for i=2:length(a.index)
            result *= ","*string(a.index[i]);
        end
        result *= "]";
    elseif a.type==IRtypes.matrix_data && length(a.index) > 1
        result *= "{Matrix("*string(a.index[1]);
        result *= ","*string(a.index[2]) * ")}";
    elseif a.type==IRtypes.vector_data && length(a.index) > 0
        result *= "{Vector("*string(a.index[1]) * ")}";
    end
    return result;
end
function IR_string(a::IR_data_access)
    result = a.type==31 ? "read(" : "write(";
    result *= string(a.data)*")";
    return result;
end

# Prints the IR_part as a tree
function print_IR_tree(a::IR_part)
    AbstractTrees.print_tree(a);
end

#############################################################################################################
include("IR_utils.jl");
include("IR_build_FEM.jl");

#############################################################################################################
# stuff to be removed

struct tmpsym
    name::String
    index::Int
end

function testIR()
    IRtypes = IR_entry_types();
    
    xalloc = IR_statement_node(IRtypes.allocate_statement, :x, IR_allocate_node(IRtypes.float_64_data, [100]));
    
    xdata = IR_data_node(IRtypes.array_data, :x, [:eid]);
    adata = IR_data_node(IRtypes.array_data, :a, [:eid]);
    bdata = IR_data_node(IRtypes.scalar_data, :b);
    
    pluseval = IR_eval_node(IRtypes.math_eval, :add, [adata,bdata]);
    assign = IR_statement_node(IRtypes.assign_statement, xdata, pluseval);
    assign2 = IR_statement_node(IRtypes.assign_statement, xdata, 1);
    equal = IR_eval_node(IRtypes.func_eval, :(==), [xdata, 0]);
    
    cond = IR_conditional_node(equal, [assign2], nothing);
    
    eloop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, 100, [assign, cond]);
    #########################################################
    xidata = IR_data_node(IRtypes.array_data, :x, [:i]);
    
    pluseval2 = IR_eval_node(IRtypes.math_eval, :add, [xidata, 2]);
    assign2 = IR_statement_node(IRtypes.assign_statement, xidata, pluseval2);
    
    iloop = IR_loop_node(IRtypes.index_loop, :i, :i, 1, 100, [assign2]);
    #########################################################
    tloop = IR_loop_node(IRtypes.time_loop, :time, :ti, 1, 75, [eloop, iloop]);
    
    display(xalloc);
    display(tloop);
    
    # display(eloop.deps);
    println("Dependencies:");
    display(tloop.deps);
    
    ###############################
    # For testing
    
    dimension = 2;
    needed_derivative_matrices = [true, true, false, false];
    geo_factors_index = :row;
    derivative_matrix_loop_body = Vector{IR_part}(undef,0);
    
    row_col_matrix_index = IR_eval_node(IRtypes.special_eval, :ROWCOL_TO_INDEX, [:row, :col, :nnodes]);
    push!(derivative_matrix_loop_body, IR_statement_node(IRtypes.assign_statement, :idx, row_col_matrix_index));
    gf_rx = IR_data_node(IRtypes.array_data, :geo_factors_rx, [:eid, geo_factors_index]);
    gf_ry = IR_data_node(IRtypes.array_data, :geo_factors_ry, [:eid, geo_factors_index]);
    gf_rz = IR_data_node(IRtypes.array_data, :geo_factors_rz, [:eid, geo_factors_index]);
    gf_sx = IR_data_node(IRtypes.array_data, :geo_factors_sx, [:eid, geo_factors_index]);
    gf_sy = IR_data_node(IRtypes.array_data, :geo_factors_sy, [:eid, geo_factors_index]);
    gf_sz = IR_data_node(IRtypes.array_data, :geo_factors_sz, [:eid, geo_factors_index]);
    gf_tx = IR_data_node(IRtypes.array_data, :geo_factors_tx, [:eid, geo_factors_index]);
    gf_ty = IR_data_node(IRtypes.array_data, :geo_factors_ty, [:eid, geo_factors_index]);
    gf_tz = IR_data_node(IRtypes.array_data, :geo_factors_tz, [:eid, geo_factors_index]);
    Qr = IR_data_node(IRtypes.array_data, :refel_Qr, [:idx]);
    Qs = IR_data_node(IRtypes.array_data, :refel_Qs, [:idx]);
    Qt = IR_data_node(IRtypes.array_data, :refel_Qt, [:idx]);
    for i=1:3
        if needed_derivative_matrices[i]
            if i==1
                gf_r = gf_rx;
                gf_s = gf_sx;
                gf_t = gf_tx;
            elseif i==2
                gf_r = gf_ry;
                gf_s = gf_sy;
                gf_t = gf_ty;
            else
                gf_r = gf_rz;
                gf_s = gf_sz;
                gf_t = gf_tz;
            end
            RQn = IR_data_node(IRtypes.array_data, Symbol("RQ"*string(i)), [:idx]);
            
            tmpexprargs = [IR_eval_node(IRtypes.math_eval, :(*), [gf_r, Qr])];
            if dimension > 1 push!(tmpexprargs, IR_eval_node(IRtypes.math_eval, :(*), [gf_s, Qs])) end
            if dimension > 2 push!(tmpexprargs, IR_eval_node(IRtypes.math_eval, :(*), [gf_t, Qt])) end
            
            tmpexpr = IR_eval_node(IRtypes.math_eval, :(+), tmpexprargs);
            
            tmpst = IR_statement_node(IRtypes.assign_statement, RQn, tmpexpr);
            push!(derivative_matrix_loop_body, tmpst);
        end
    end
    loop_col_nodes = IR_loop_node(IRtypes.space_loop, :local_nodes, :col, 0, 0, derivative_matrix_loop_body);
    loop_row_qnodes = IR_loop_node(IRtypes.space_loop, :quad_nodes, :row, 0, 0, [loop_col_nodes])
    
    display(loop_row_qnodes)
    
    ##################################
    
    derivative_matrix_building = Vector{IR_part}(undef,0);
    for dim=1:3
        if needed_derivative_matrices[dim]
            push!(derivative_matrix_building, 
                  IR_statement_node(IRtypes.call_statement, nothing, IR_eval_node(IRtypes.func_eval, :build_derivative_matrix, [dim, 1, Symbol("RQ"*string(dim))])));
        end
    end
    
    display(derivative_matrix_building)
    
    #######################################
    
    row_col_matrix_index = IR_eval_node(IRtypes.special_eval, :ROWCOL_TO_INDEX, [:row, :col, :nnodes]);
    array_zero_loop = IR_loop_node(IRtypes.space_loop, :nodes, :row, 0, 0, [
            IR_loop_node(IRtypes.space_loop, :nodes, :row, 0, 0, [
                IR_statement_node(IRtypes.assign_statement, 
                    IR_data_node(IRtypes.array_data, :elemental_matrix, [row_col_matrix_index]),
                    0.0)
            ]),
            IR_statement_node(IRtypes.assign_statement, 
                IR_data_node(IRtypes.array_data, :elemental_vector, [:row]),
                0.0)
        ]);
    
    display(array_zero_loop)
    
    ############################################
    
    # allocate coefficient vectors
    allocate_part = Vector{IR_part}(undef,0);
    needed_coefficient_vectors = [tmpsym("f", 1)];
    for c in needed_coefficient_vectors
        nodal_coef_name = "NODALvalue_" * c.name * "_" * string(c.index);
        quad_coef_name = "value_" * c.name * "_" * string(c.index);
        
        push!(allocate_part, IR_statement_node(IRtypes.allocate_statement,
            IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name)),
            IR_allocate_node(IRtypes.float_64_data, [:nodes_per_element])));
        push!(allocate_part, IR_statement_node(IRtypes.allocate_statement,
            IR_data_node(IRtypes.array_data, Symbol(quad_coef_name)),
            IR_allocate_node(IRtypes.float_64_data, [:quadnodes_per_element])));
    end
    
    # Compute all coefficients
    # First get x,y,z,t,nodeID
    coef_loop_body = [
        # nodeID = mesh_loc2glb[eid, ni]
        IR_statement_node(IRtypes.assign_statement, :nodeID, 
            IR_data_node(IRtypes.array_data, :mesh_loc2glb, [:eid, :ni])),
        # t = 0.0
        IR_statement_node(IRtypes.assign_statement, :t, 0.0),
        # x = mesh_allnodes[nodeID*dimension]
        IR_statement_node(IRtypes.assign_statement, :x, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension])
            ]))
    ];
    if dimension > 1
        # y = mesh_allnodes[nodeID*dimension+1]
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :y, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(+), [IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension]), 1])
            ])))
    end
    if dimension > 2
        # z = mesh_allnodes[nodeID*dimension+2]
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :z, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(+), [IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension]), 2])
            ])))
    end
    
    # Evaluate coefficient functions
    for c in needed_coefficient_vectors
        nodal_coef_name = "NODALvalue_" * c.name * "_" * string(c.index);
        genfunction_index = 2; # Need to determine this...
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement,
            IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:ni]),
            IR_eval_node(IRtypes.special_eval, :GENFUNCTION, [genfunction_index, :x, :y, :z, :t, :nodeID])));
    end
    
    # Build the coef eval loop
    coef_eval_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 0, 0, coef_loop_body);
    
    # Interpolate, change to modal, apply derivatives if needed
    coef_init_loop_body = [];
    coef_interp_loop_body = [];
    row_col_matrix_index = IR_eval_node(IRtypes.special_eval, :ROWCOL_TO_INDEX, [:row, :col, :nnodes]);
    for c in needed_coefficient_vectors
        quad_coef_name = "value_" * c.name * "_" * string(c.index);
        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(quad_coef_name), [:row]);
        nodal_coef_name = "NODALvalue_" * c.name * "_" * string(c.index);
        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:col]);
        refelQ = IR_data_node(IRtypes.array_data, :refel_Q, [row_col_matrix_index]);
        push!(coef_init_loop_body, IR_statement_node(IRtypes.assign_statement, quad_coef_node, 0.0));
        push!(coef_interp_loop_body, IR_statement_node(IRtypes.assign_statement,
            quad_coef_node,
            IR_eval_node(IRtypes.math_eval, :(+), [quad_coef_node, IR_eval_node(IRtypes.math_eval, :(*), [refelQ, nodal_coef_node])])));
    end
    
    push!(coef_init_loop_body, IR_loop_node(IRtypes.space_loop, :nodes, :row, 0, 0, coef_interp_loop_body))
    
    coef_interp_loop = IR_loop_node(IRtypes.space_loop, :qnodes, :row, 0, 0, coef_init_loop_body);
    
    display(coef_eval_loop)
    display(coef_interp_loop);
    
    everythingvec::Vector{IR_part} = [tloop, loop_row_qnodes, array_zero_loop, coef_eval_loop, coef_interp_loop];
    println(typeof(everythingvec))
    everything = IR_block_node(everythingvec);
    print_tree(everything, maxdepth=20)
    
    # pyplot()
    # plot(TreePlot(coef_eval_loop), method=:tree, fontsize=10, nodeshape=:rect)
    
    
end

end # module


# using .IntermediateRepresentation

# testIR()