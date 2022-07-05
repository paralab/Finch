#=
This is a macro used by modules to import various Finch items.
=#

finch_symbols_to_import = [
    #Constants
    :JULIA, :CPP, :MATLAB, :DENDRO, :HOMG, :CUSTOM_GEN_TARGET,
    :SQUARE, :IRREGULAR, :UNIFORM_GRID, :TREE, :UNSTRUCTURED, 
    :CG, :DG, :HDG, :FV, :MIXED,
    :NODAL, :MODAL, :CELL, :LEGENDRE, :UNIFORM, :GAUSS, :LOBATTO, 
    :NONLINEAR_NEWTON, :NONLINEAR_SOMETHING, 
    :EULER_EXPLICIT, :EULER_IMPLICIT, :CRANK_NICHOLSON, :RK4, :LSRK4, :PECE, :MIXED_STEPPER,
    :DEFAULT_SOLVER, :PETSC_SOLVER, :CUDA_SOLVER,
    :VTK, :RAW_OUTPUT, :CUSTOM_OUTPUT, 
    :DIRICHLET, :NEUMANN, :ROBIN, :NO_BC, :FLUX, :SYMMETRIC,
    :MSH_V2, :MSH_V4,
    :SCALAR, :VECTOR, :TENSOR, :SYM_TENSOR, :VAR_ARRAY,
    :LHS, :RHS,
    :LINEMESH, :QUADMESH, :TRIMESH, :HEXMESH,
    
    #Structs
    :Variable, :Coefficient, :GenFunction, :CallbackFunction, :Indexer, :Grid, :Refel, :Finch_prob, :GeometricFactors, :VariableTransform,
    :SymExpression, :SymEntity, :CacheSim,
    
    #Functions
    :geometric_factors, :geometric_factors_face, :build_deriv_matrix, :evaluate_coefficient, :exchange_ghosts,
    
    #Log and profiling
    :log_entry, :printerr, :timer_output, :timeit,
    
    #Data
    :config, :prob, :variables, :coefficients, :parameters, :test_functions, :indexers, :variable_transforms,
    :mesh_data, :grid_data, :refel, :time_stepper, :language, :gen_framework,
    :elemental_order, :genfunctions, :callback_functions, :geo_factors, :use_cachesim, :custom_gen_funcs,
    
    #Specific to FV
    :FVInfo, :fv_info, :FV_cell_to_node, :FV_node_to_cell, :FV_reconstruct_value_left_right,
    :ParentMaps, :parent_maps, :fv_grid, :fv_refel, :fv_geo_factors, :build_local_patch,
    
    # Packages? should this be done here?
    :MPI, :PETSc, :TimerOutputs
];

macro import_finch_symbols()
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