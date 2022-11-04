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
- comment - a string that is simply a comment
=#

module IntermediateRepresentation

# Implement AbstractTrees for plotting/printing/using other tools if desired
# using AbstractTrees
try
    using AbstractTrees
catch e
    println("AbstractTrees package is not yet installed. Installing now.");
    using Pkg
    Pkg.add("AbstractTrees")
    using AbstractTrees
end

export IR_entry_types, IR_string, print_IR, repr_IR
export IR_part, IR_data_node, IR_data_access, IR_operation_node, IR_block_node, IR_loop_node, IR_conditional_node, IR_comment_node
export build_IR_fem, build_IR_fvm

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

AbstractTrees.children(a::Nothing) = ();
AbstractTrees.children(a::Symbol) = ();
# AbstractTrees.children(a::Number) = ();
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
    index_loop::Int8  # = 4 # other indices in a for loop
    while_loop::Int8  # = 5 # while loop
    
    # operation types
    allocate_op::Int8 # = 11 # allocation has at least two args: type, dims...
    assign_op::Int8     # = 12 # assignment is like lhs = rhs, lhs is a data node, rhs is data or op
    function_op::Int8   # = 13 # function evaluation has name and tuple of args like name(args)
    math_op::Int8       # = 14 # a function_op for arithmetic just in case it is useful to specify
    member_op::Int8     # = 15 # access a member of a struct
    named_op::Int8      # = 16 # a special op keyword that will be interpreted by codegen as needed
    
    # data types
    int_data::Int8      # = 21 # These are for int and float of unspecified size
    float_data::Int8    # = 22 #
    
    int_32_data::Int8    # = 25
    int_64_data::Int8    # = 26
    float_32_data::Int8  # = 27
    float_64_data::Int8  # = 28
    boolean_data::Int8   # = 29
    
    # access types
    read_access::Int8  # = 31 
    write_access::Int8 # = 32
    
    IR_entry_types() = new(
        Dict{Int8, String}(
            1=>"time loop",
            2=>"space loop",
            3=>"dof loop",
            4=>"index loop",
            5=>"while loop",
            
            11=>"allocate op",
            12=>"assign op",
            13=>"function op",
            14=>"math op",
            15=>"member op",
            16=>"named op",
            
            21=>"Int",
            22=>"Float64",
            25=>"Int32",
            26=>"Int64",
            27=>"Float32",
            28=>"Float64",
            29=>"Bool",
            
            31=>"read",
            32=>"write"
        ),
        1,2,3,4,5, 11,12,13,14,15,16, 21,22,25,26,27,28,29, 31,32
        );
end

# This represents a location for data.
# It could be a variable name or an array name with index symbols.
struct IR_data_node <: IR_part
    type::Int8 # The type of data
    label::Symbol
    size::Vector{Union{Symbol,Int,IR_part}} # scalars = [], unknown array size = [:?]
    index::Vector{Union{Symbol,Int,IR_part}} # only used when refering to an element of an array
    IR_data_node(t::Int8, v::Symbol) = new(t,v,[],[]);
    IR_data_node(t::Int8, v::Symbol, s::Vector, i::Vector) = new(t,v,s,i);
end
AbstractTrees.children(a::IR_data_node) = a.index;
AbstractTrees.printnode(io::IO, a::IR_data_node) = print(io, string(a.label));

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
    elseif a.type == IRtypes.math_op
        return a.args[2:end];
    elseif a.type == IRtypes.member_op
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
        print(io, string(a.args[1]));
    elseif a.type == IRtypes.math_op
        print(io, string(a.args[1]));
    elseif a.type == IRtypes.member_op
        print(io, string(a.args[1]));
    elseif a.type == IRtypes.named_op
        print(io, string(a.args[1]));
    else
        print(io, "unknown op");
    end
end

# This is a block of IR_parts such as a set of statements.
# Use for the body of a loop or conditional or for a top-level container
mutable struct IR_block_node <: IR_part
    parts::Vector{IR_part}
    label::String
    IR_block_node(p::Vector) = new(p, "")
    IR_block_node(p::Vector, s::String) = new(p, s)
end
AbstractTrees.children(a::IR_block_node) = a.parts;
AbstractTrees.printnode(io::IO, a::IR_block_node) = print(io,"block:"*a.label);

# This is a loop.
# The type could be time, space, dof, or index or while
# It has a name of the collection it is looping over, the iterator symbol, and limits.
# If limits are both zero it loops over the whole collection.
# If it is a while loop, the "last" component is an IR_part for a condition,
# and the iterator is included and incremented at the beginning of the loop.
mutable struct IR_loop_node <: IR_part
    type::Int8   # The type of loop
    collection::Symbol
    iterator::Symbol
    first::Union{Int,Symbol,IR_part}
    last::Union{Int,Symbol,IR_part}
    body::IR_block_node #Vector{IR_part}
    deps::Vector{IR_data_access}
    function IR_loop_node(t::Int8, c::Symbol, i::Symbol, f, l, b)
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
    elseif a.type == IRtypes.while_loop
        str *= "while "*string(a.last);
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
    function IR_conditional_node(c, b)
        IR_conditional_node(c, b, nothing)
    end
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

# This is just a comment string
struct IR_comment_node <: IR_part
    string::String # the comment
end
AbstractTrees.children(a::IR_comment_node) = ();
AbstractTrees.printnode(io::IO, a::IR_comment_node) = print(io, "(comment)"*a.string);

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
    if !(a.elsepart === nothing) && length(a.elsepart.parts) > 0
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
        
    elseif a.type == IRtypes.member_op
        result = string(a.args[1])*"."*string(a.args[2]);
        
    else
        result = "";
    end
    
    return result;
end
function IR_string(a::IR_data_node)
    IRtypes = IR_entry_types();
    result = string(a.label);
    if length(a.size)>0 && length(a.index) > 0
        result *= "["*string(a.index[1]);
        for i=2:length(a.index)
            result *= ","*string(a.index[i]);
        end
        result *= "]";
    end
    return result;
end
function IR_string(a::IR_data_access)
    result = a.type==31 ? "read(" : "write(";
    result *= string(a.data)*")";
    return result;
end
function IR_string(a::IR_comment_node)
    return "# "*a.string;
end

# Prints the IR_part as a tree
function print_IR(a::IR_part)
    AbstractTrees.print_tree(a);
end
function repr_IR(a::IR_part)
    return AbstractTrees.repr_tree(a, maxdepth=20);
end

#############################################################################################################
include("IR_utils.jl");
include("IR_build_shared.jl");
include("IR_build_FEM.jl");
include("IR_build_FVM.jl");

#############################################################################################################

end # module
