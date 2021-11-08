
# External modules
using SparseArrays
using LinearAlgebra
using Random

######################
# NOTE: This is not a long term solution.
# Once the package is set up, we can put this dependency in the .toml
try
    using SymEngine
catch e
    println("Julia SymEngine is not yet installed. Installing now.");
    using Pkg
    Pkg.add("SymEngine")
    using SymEngine
end
try
    using WriteVTK
catch e
    println("WriteVTK is not yet installed. Installing now.");
    using Pkg
    Pkg.add("WriteVTK")
    using WriteVTK
end
# try
#     using Latexify
# catch e
#     println("Latexify is not yet installed. Installing now.");
#     using Pkg
#     Pkg.add("Latexify")
#     using Latexify
# end
######################

# include these first
include("finch_constants.jl");
include("finch_config.jl");
include("finch_prob.jl");
include("finch_import_symbols.jl");

# include these next
include("logging.jl");

include("refel.jl");
include("mesh_data.jl");
include("mesh_read.jl");
include("mesh_write.jl")
include("grid.jl");
include("simple_mesh.jl");
include("recursive_ordering.jl");
include("tiled_ordering.jl");
include("ef_ordering.jl");
include("bad_ordering.jl");

include("general_utils.jl");
include("function_utils.jl");
include("symexpression.jl");
include("variables.jl");
include("coefficient.jl");
include("parameter.jl");
include("indexer.jl");
include("geometric_factors.jl");
include("dg_utils.jl");
include("time_steppers.jl");
include("fv_utils.jl");
include("fv_neighborhood.jl");
include("grid_parent_child.jl");
include("polyharmonic_interp.jl");

include("output_data.jl");
include("cachsim_solve.jl");

# Finch submodules
include("cachesim/CachesimOut.jl");
using .CachesimOut
include("SymbolicParser.jl")
using .SymbolicParser
include("CodeGenerator.jl");
using .CodeGenerator
include("DGSolver.jl");
using .DGSolver
include("CGSolver.jl");
using .CGSolver
include("FVSolver.jl");
using .FVSolver
