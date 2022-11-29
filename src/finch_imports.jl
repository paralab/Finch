
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
            exit(0);
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
        println("You will need to restart Julia.");
        try
            Pkg.add("MPIPreferences")
            using MPIPreferences
            MPIPreferences.use_system_binary()
            global need_restart = true;
        catch e2
            println("There was a problem setting up MPI.jl. Please do this manually before attempting a parallel run.");
            println("Continuing without MPI.");
        end
        
    else
        println("MPI package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("MPI")
            using MPI
            println("Note that MPI may require some manual configuration to work properly.")
            println("Now attempting to build the MPI package using an MPI implementation installed on your system.");
            println("You will need to restart Julia.");
            try
                Pkg.add("MPIPreferences")
                using MPIPreferences
                MPIPreferences.use_system_binary()
                global need_restart = true;
            catch e2
                println("There was a problem setting up MPI.jl. Please do this manually before attempting a parallel run.");
                println("Continuing without MPI.");
            end
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

try
    using LinearMaps
catch e
    if force_package_install
        println("LinearMaps package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("LinearMaps")
        using LinearMaps
        
    else
        println("LinearMaps package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("LinearMaps")
            using LinearMaps
        else
            println("Continuing without LinearMaps. Note that matrix-free solvers will not be available.")
        end
    end
end

try
    using IncompleteLU
catch e
    if force_package_install
        println("IncompleteLU package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("IncompleteLU")
        using IncompleteLU
        
    else
        println("IncompleteLU package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("IncompleteLU")
            using IncompleteLU
        else
            println("Continuing without IncompleteLU.")
        end
    end
end

try
    using AlgebraicMultigrid
catch e
    if force_package_install
        println("AlgebraicMultigrid package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("AlgebraicMultigrid")
        using AlgebraicMultigrid
        
    else
        println("AlgebraicMultigrid package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("AlgebraicMultigrid")
            using AlgebraicMultigrid
        else
            println("Continuing without AlgebraicMultigrid.")
        end
    end
end

try
    using IterativeSolvers
catch e
    if force_package_install
        println("IterativeSolvers package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("IterativeSolvers")
        using IterativeSolvers
        
    else
        println("IterativeSolvers package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("IterativeSolvers")
            using WriteVTK
        else
            println("Continuing without IterativeSolvers. Note that iterative solvers will not be available.")
        end
    end
    
end

try
    using PETSc
catch e
    # println("PETSc package is not yet installed. It is optional and may cause issues depending on your system.\nWould you like to install now? y/n.");
    # response = readline();
    # if response=="Y" || response=="y"
    #     using Pkg
    #     Pkg.add("PETSc")
    #     try
    #         using PETSc
    #     catch e
    #         println("There was an issue with using PETSc. Check the Julia PETSc.jl documentation to install correctly.")
    #         exit(0);
    #     end
    # else
    #     println("Continuing without PETSc.")
    #     PETSc = nothing;
    # end
    println("PETSc package is not available. It is optional, but PETSc solvers can't be used.")
    
end

try
    using TimerOutputs
catch e
    if force_package_install
        println("TimerOutputs package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("TimerOutputs")
        using TimerOutputs
        
    else
        println("TimerOutputs package is not yet installed. It is required.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("TimerOutputs")
            using TimerOutputs
        else
            println("It is a required package. Exiting Finch.")
            exit(0);
        end
    end
end

try
    using Zygote
catch e
    if force_package_install
        println("Zygote package is not yet installed. Installing now.");
        using Pkg
        Pkg.add("Zygote")
        using Zygote
        
    else
        println("Zygote package is not yet installed. It is optional.\nWould you like to install now? y/n.");
        response = readline();
        if response=="Y" || response=="y"
            using Pkg
            Pkg.add("Zygote")
            using Zygote
        else
            println("Continuing without Zygote")
        end
    end
end

# Implement AbstractTrees for plotting/printing/using other tools if desired
try
    using AbstractTrees
catch e
    println("AbstractTrees package is not yet installed. Installing now.");
    using Pkg
    Pkg.add("AbstractTrees")
    using AbstractTrees
end
