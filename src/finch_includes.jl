#=
Finch code files to include.
=#

# include these first
# include("finch_config.jl");
# include("finch_problem.jl");
include("finch_import_symbols.jl");

# include these next
include("logging.jl");

include("refel.jl");
include("mesh_data.jl");
include("mesh_read.jl");
include("mesh_write.jl")
include("grid.jl");
include("grid_parent_child.jl");
include("subdomain.jl");
include("simple_mesh.jl");
include("recursive_ordering.jl");
include("tiled_ordering.jl");
include("ef_ordering.jl");
include("bad_ordering.jl");
include("partition.jl");

include("general_utils.jl");
include("function_utils.jl");
include("symexpression.jl");
include("basic_ops.jl");
include("variables.jl");
include("coefficient.jl");
# include("parameter.jl");
# include("indexer.jl");
include("geometric_factors.jl");
include("dg_utils.jl");
include("time_steppers.jl");
include("fv_utils.jl");
include("fv_neighborhood.jl");

include("polyharmonic_interp.jl");
include("solver_utils.jl");
include("boundary_utils.jl");
include("parallel_utils.jl");

include("IntermediateRepresentation.jl")

include("output_data.jl");

include("macros.jl");
include("finch_interface.jl");

# cachesim parts
include("cachesim/CacheSim.jl");
using .CacheSim
include("cachesim_utils.jl");
include("cachesim_solve.jl");

# Finch submodules
include("SymbolicParser.jl")
using .SymbolicParser

include("CodeGenerator.jl");
using .CodeGenerator
