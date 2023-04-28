#=
This is a macro used by modules to import various Finch symbols.
=#
macro import_finch_symbols()
    finch_symbols_to_import = [
        #Constants
        :JULIA, :CPP, :MATLAB, :FINCH, :DENDRO, :HOMG, :CUSTOM_GEN_TARGET,
        :SQUARE, :IRREGULAR, :UNIFORM_GRID, :TREE, :UNSTRUCTURED, 
        :CG, :DG, :HDG, :FV, :MIXED,
        :NODAL, :MODAL, :CELL, :LEGENDRE, :UNIFORM, :GAUSS, :LOBATTO, 
        :NONLINEAR_NEWTON, :NONLINEAR_SOMETHING, 
        :EULER_EXPLICIT, :EULER_IMPLICIT, :CRANK_NICHOLSON, :RK4, :LSRK4, :BDF2, :PECE, :MIXED_STEPPER,
        :DEFAULT_SOLVER, :PETSC_SOLVER, :CUDA_SOLVER,
        :VTK, :RAW_OUTPUT, :CUSTOM_OUTPUT, 
        :DIRICHLET, :NEUMANN, :ROBIN, :NO_BC, :FLUX, :SYMMETRIC,
        :MSH_V2, :MSH_V4, :METIS, :FENNEL,
        :SCALAR, :VECTOR, :TENSOR, :SYM_TENSOR, :VAR_ARRAY,
        :LHS, :RHS,
        :LINEMESH, :QUADMESH, :TRIMESH, :HEXMESH,
        
        #Structs
        :FinchState, :Variable, :Coefficient, :GenFunction, :CallbackFunction, :Indexer, :Grid, 
        :Refel, :FinchProblem, :GeometricFactors, :VariableTransform, :Subdomain,
        :SymExpression, :SymEntity, :SymOperator, :CacheSim, :FVInfo, :ParentMaps,
        
        #Functions
        :geometric_factors, :evaluate_coefficient,
        
        #Log and profiling
        :log_entry, :printerr, :timeit,
        
        #Data
        :finch_state,
        
        #Specific to FV
        :FV_cell_to_node, :FV_node_to_cell, :FV_reconstruct_value_left_right,
        :build_local_patch,
        
        # Packages? should this be done here?
        :LinearMaps, :IterativeSolvers, :AlgebraicMultigrid, :IncompleteLU, :MPI, :TimerOutputs, :Zygote
    ];
    
    
    # Form an Expr that looks like: 
    #     import ..Finch: a, b, c, d, ...
    num_symbols = length(finch_symbols_to_import);
    ex = :(import ..Finch: a);
    symb = Array{Expr,1}(undef,num_symbols);
    for i=1:num_symbols
        symb[i] = Expr(:., finch_symbols_to_import[i]);
    end
    ex.args[1].args = [ex.args[1].args[1]];
    append!(ex.args[1].args, symb);
    
    ex
end