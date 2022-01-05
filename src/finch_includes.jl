
# External modules
using SparseArrays
using LinearAlgebra
using Random

######################
# NOTE: This is not a long term solution.
# Once the package is set up, we can put this dependency in the .toml
force_package_install = true;
need_restart = false;
try
    using SymEngine
catch e
    if force_package_install
        println("SymEngine package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("SymEngine")
        using SymEngine
        
    else
        println("SymEngine package is not yet installed. It is required.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("SymEngine")
            using SymEngine
        else
            println("It is a required package. Exiting Finch.")
        end
    end
    
end
try
    using WriteVTK
catch e
    if force_package_install
        println("WriteVTK package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("WriteVTK")
        using WriteVTK
        
    else
        println("WriteVTK package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("WriteVTK")
            using WriteVTK
        else
            println("Continuing without WriteVTK. Note that VTK output will not be possible.")
        end
    end
    
end
try
    using MPI
catch e
    if force_package_install
        println("MPI package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("MPI")
        using MPI
        println("Note that MPI may require some manual configuration to work properly.")
        println("Now attempting to build the MPI package using an MPI implementation installed on your system.");
        ENV["JULIA_MPI_BINARY"] = "system";
        Pkg.build("MPI"; verbose=true);
        need_restart = true;
        
    else
        println("MPI package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("MPI")
            using MPI
            println("Note that MPI may require some manual configuration to work properly.")
            println("Now attempting to build the MPI package using an MPI implementation installed on your system.");
            ENV["JULIA_MPI_BINARY"] = "system";
            Pkg.build("MPI"; verbose=true);
            need_restart = true;
        else
            println("Continuing without MPI.")
        end
    end
    
end
if @isdefined(MPI)
    try
        using METIS_jll
    catch e
        if force_package_install
            println("METIS_jll package is not yet installed. Installing now.");
            using Pkg
            Pkg.add("METIS_jll")
            using METIS_jll: libmetis
            
        else
            println("METIS_jll package is not yet installed. It is optional, but needed for effective MPI parallelization.\nWould you like to install now? y/n.");
            response = readline();
            if response=="Y" || response=="y"
                using Pkg
                Pkg.add("METIS_jll")
                using METIS_jll: libmetis
            else
                println("Continuing without METIS_jll.")
            end
        end
        
    end
end

if need_restart
    println("Julia may need to restart to use newly installed packages. Exiting now. Please run again.")
    exit(0);
end
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
include("partition.jl");
include("ghost_exchange.jl");

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
