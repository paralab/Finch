#=
This contains all of the pieces needed to add a new code gen target.
The following three functions must be provided.

1. get_external_language_elements() - file extensions, comment chars etc.
2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code

You will have access to anything in the Finch module scope because this file will be included there.
=#

function get_external_language_elements()
    file_extension = ".cpp";
    comment_char = "//";
    block_comment = ["/*"; "*/"];
    
    return (file_extension, comment_char, block_comment);
end

#=
Translate the symbolic layer into code for the elemental assembly.
var - unknown variable or array of variables.
entities - leaf nodes of the computational graph. Type is SymEntity(see symexpression.jl) or number.
terms - array of additive terms, separated for convenience. [t1, t2, t3] for expression t1+t2+t3 (see note below)
lorr - LHS or RHS will have constant values LHS or RHS which are "lhs", "rhs"
vors - volume or surface will have values "volume" or "surface"

Note on terms: Each term is a SymExpression(see symexpression.jl) which is essentially a Julia Expr and list of entities.
The Expr holds the computational graph for each term. LHS and RHS are determined by separating terms that involve 
unknowns from those that don't. Also, in well formed FEM problems each term should be multiplied by a test function,
and thus quadrature can be done for each term separately.
=#
function generate_external_code_layer(var, entities, terms, lorr, vors)
    # This can be divided up into smaller functions however you wish.
    # The whole body of code is stored in a string that will be inserted in the 
    # appropriate function during file writing.
    code = "";
    
    # Allocate/Evaluate or fetch the values for each needed entity.
    # Note: This could be staged and interleaved with the calculation part for complex problems. TODO
    (prep_code, to_delete) = dendrotarget_prepare_needed_values(entities, var, lorr, vors);
    code *= prep_code;
    code *= "\n";
    
    # Form the final elemental calculation
    code *= dendrotarget_make_elemental_computation(terms, to_delete, var, lorr, vors);
    
    return code;
end

# Create all of the code files.
# The input are code strings generated above. They will be inserted into their appropriate files.
function generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf)
    if codegen_params === nothing || length(codegen_params) < 5
        # default values (max_depth=6, wavelet_tol = 0.1, partition_tol = 0.3, solve_tol = 1e-6, max_iters = 100)
        params=(5, 1, 0.3, 0.000001, 100);
    else
        params = codegen_params;
    end
    # Again, this can be split up into smaller functions as you wish.
    dendro_utils_file();            # various math and other utils that will be filled as needed.
    dendro_main_file(var);          # The main file
    dendro_config_file(params);     # Dendro parameters and other configuration
    dendro_prob_file();             # Problem details: boundary conditions, initial conditions, etc.
    dendro_genfunction_file();      # Functions for coefficients, boundary conditions etc.
    dendro_bilinear_file(lhs_vol);  # The matrix part (feMatrix)
    dendro_linear_file(rhs_vol);    # The vector part (feVector)
    if Finch.prob.time_dependent
        dendro_stepper_file();      # time stepper
    end
    dendro_output_file();           # Output code
    dendro_cmake_file();            # Writes a cmake file
    dendro_readme_file();           # readme
end

# 
global octree_building_function = nothing;
function build_octree_with(a)
    if typeof(a) == Coefficient && typeof(a.value[1]) == GenFunction
        global octree_building_function = a.value[1];
    elseif typeof(a) == GenFunction
        global octree_building_function = a;
    elseif typeof(a) == String
        for c in coefficients
            if a == string(c.symbol) && typeof(c.value[1]) == GenFunction
                global octree_building_function = c.value[1];
            end
        end
        if octree_building_function === nothing
            printerr("Unexpected input to build_octree_with(). Using default.")
        end
    else
        printerr("Unexpected input to build_octree_with(). Using default.")
    end
end
export build_octree_with;

#########################################################
# code writing utilities
#########################################################

# numbers: 2 -> "2"
# strings: "thing" -> "\"thing\""
# arrays: [1 2; 3 4] -> "{{1,2},{3,4}}"
function cpp_gen_string(v)
    if typeof(v) == String
        return "\""*v*"\"";
        
    elseif typeof(v) <: Number
        return string(v);
        
    elseif typeof(v) <: Array
        if ndims(v) == 1 # "{a, b, c}"
            n = length(v);
            str = "{";
            for i=1:n
                str = str*cpp_gen_string(v[i]);
                if i < n
                    str = str*", ";
                end
            end
            str = str*"}";
        elseif ndims(v) == 2 # "{{a, b}, {c, d}}
            (n,m) = size(v);
            str = "{";
            for i=1:n
                for j=1:m
                    str = str*cpp_gen_string(v[i,j]);
                    if i*j < n*m
                        str = str*", ";
                    end
                end
            end
            str = str*"}";
        end
        return str;
        
    elseif typeof(v) == GenFunction
        return v.name;
    else
        return string(v);
    end
end

#= Generates:
std::function<returntype(argtypes)> name = [captures](args){
    (content)
};
# Note: content should be an array of lines to allow indentation
=#
function cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content)
    lines = fill("", length(content)+2);
    inner_indent = indent*"    ";
    
    arg = string(argtypes[1])*" "*string(args[1]);
    argtype = string(argtypes[1]);
    for i=2:length(args)
        arg = arg*", "*string(argtypes[i])*" "*string(args[i]);
        argtype = argtype*", "*string(argtypes[i]);
    end
    ret = string(ret);
    
    lines[1] = indent*"std::function<"*rettype*"("*argtype*")> "*name*" = ["*captures*"]("*arg*"){";
    for i=1:length(content)
        lines[i+1] = inner_indent*content[i];
    end
    lines[end] = indent*"};";
    
    return lines;
end

#= Generates:
returntype name(args){
    (content)
}
# Note: content should be an array of lines to allow indentation
=#
function cpp_function_def(indent, name, args, argtypes, rettype, content)
    lines = fill("", length(content)+2);
    inner_indent = indent*"    ";
    
    arg = string(argtypes[1])*" "*string(args[1]);
    argtype = string(argtypes[1]);
    for i=2:length(args)
        arg = arg*", "*string(argtypes[i])*" "*string(args[i]);
        argtype = argtype*", "*string(argtypes[i]);
    end
    
    lines[1] = indent*rettype*" "*name*"("*arg*"){";
    for i=1:length(content)
        lines[i+1] = inner_indent*content[i];
    end
    lines[end] = indent*"}";
    
    return lines;
end

#= Generates:
for(iterator=range[1]; iterator<range[2]; iterator+=step){
    (content)
}
# Note: content should be an array of lines to allow indentation
=#
function cpp_for_loop(indent, iterator, range, step, content)
    lines = [];
    inner_indent = indent*"    ";
    push!(lines, indent*"for("*iterator*" = "*string(range[1])*";"*iterator*" < "*string(range[2])*";"*iterator*" += "*string(step)*"){");
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"}")
    
    return lines;
end

# changes symbol "a" to symbol "b" in expression ex
function dendro_swap_symbol(a, b, ex)
    if typeof(ex) == Symbol
        if ex === a
            return b;
        else
            return ex;
        end
    elseif typeof(ex) <: Number
        return ex;
    elseif typeof(ex) == Expr && length(ex.args) > 1
        swapped = copy(ex);
        for i=1:length(ex.args)
            swapped.args[i] = dendro_swap_symbol(a,b,ex.args[i]);
        end
        return swapped;
    else
        return ex;
    end
end

# changes the math operators in the expr
function dendro_change_math_ops(ex)
    if typeof(ex) == Expr
        if ex.head === :.
            # a broadcast operator, change to call
            ex.head = :call
            ex.args = [ex.args[1]; ex.args[2].args];
        end
        for i=1:length(ex.args)
            ex.args[i] = dendro_change_math_ops(ex.args[i]);
        end
    elseif typeof(ex) == Symbol
        if ex === :.+   ex = :+; end
        if ex === :.-   ex = :-; end
        if ex === :.*   ex = :*; end
        if ex === :./   ex = :/; end
        if (ex === :.^ || ex === :^)   ex = :pow; end
        if ex === :abs   ex = :fabs; end
    end
    return ex;
end

# builds a c++ functional for a constant
function dendro_number_to_function(name, val)
    indent = "";
    args = ["x"; "y"; "z"; "var"];
    argtypes = ["double"; "double"; "double"; "double*"];
    ret = "";
    rettype = "void";
    captures = "gridX_to_X,gridY_to_Y,gridZ_to_Z";
    content = ["var[0] = "*string(val)*";"];
    return cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content);
end

function dendro_genfunction_to_string(genfun; xsym=:(gridX_to_X(x)), ysym=:(gridY_to_Y(y)), zsym=:(gridZ_to_Z(z)))
    newex = dendro_swap_symbol(:x, xsym, genfun.expr); # swap x for xsym
    newex = dendro_swap_symbol(:y, ysym, newex); # swap y for ysym
    newex = dendro_swap_symbol(:z, zsym, newex); # swap z for zsym
    newex = dendro_swap_symbol(:pi, :M_PI, newex); # swap pi for M_PI
    newex = dendro_change_math_ops(newex); # change operators to match C++
    s = string(newex);
    ns = replace(s, r"([\d)])([(A-Za-z])" => s"\1*\2"); # explicitly multiply with "*"
    return ns;
end

#######################################################
# Write code files

#######  ######  ##        #######   #####
##         ##    ##        ##       ###   #
#######    ##    ##        ######     ### 
##         ##    ##        ##       #   ###
##       ######  ########  #######   #####

#######################################################

function dendro_main_file(var)
    file = add_generated_file(project_name*".cpp", dir="src");
    
    valvec = string(variables[1].symvar[1]);
    
    alloc_dof = "";
    set_dof = "";
    init_fun = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar)
            nam = string(variables[i].symvar[j]);
            # if initial condition, use it, otherwise zero
            if length(prob.initial) >= i && !(prob.initial[i] === nothing)
                if typeof(prob.initial[i]) <: Array
                    init_fun = "special_"*cpp_gen_string(prob.initial[i][j]);
                else
                    init_fun = "special_"*cpp_gen_string(prob.initial[i]);
                end
            else
                init_fun = "zero_init"
            end
            alloc_dof = alloc_dof*"double * "*nam*"=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI"*nam*", false,false);\n";
            set_dof = set_dof*"octDA->setVectorByFunction("*nam*","*init_fun*",false,false,1);\n";
            if prob.time_dependent
                alloc_dof = alloc_dof*"double * "*nam*"_next=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI"*nam*"_NEXT, false,false);\n";
                set_dof = set_dof*"octDA->setVectorByFunction("*nam*"_next,"*init_fun*",false,false,1);\n";
            end
        end
    end
    # # Now coefficients are not stored as dofs
    # for i=1:length(coefficients)
    #     for j=1:length(coefficients[i].symvar)
    #         nam = string(coefficients[i].symvar[j]);
    #         fun = cpp_gen_string(coefficients[i].value[j]);
    #         alloc_dof = alloc_dof*"double * "*nam*"=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI"*nam*", false,false);\n";
    #         set_dof = set_dof*"octDA->setVectorByFunction("*nam*","*fun*",false,false,1);\n";
    #     end
    # end
    
    # functionals used for initially setting values and refining the mesh
    functional_defs = "";
    for i = 1:length(genfunctions)
        str = dendro_genfunction_to_string(genfunctions[i]);
        functional_defs *= "
    std::function<void(double,double,double,double*)> special_"*genfunctions[i].name*" = [gridX_to_X, gridY_to_Y, gridZ_to_Z](const double x,const double y,const double z,double *var){
        var[0] = FinchDendroGenfunctions::"*genfunctions[i].name*"(x,y,z,0.0);
    };";
    end
    
    # What to base the initial mesh refinement on
    if !(octree_building_function === nothing)
        mesh_base = "special_"*cpp_gen_string(octree_building_function);
    elseif length(coefficients) > 0 && typeof(coefficients[1].value[1]) == GenFunction
        mesh_base = "special_"*cpp_gen_string(coefficients[1].value[1]);
    elseif prob.time_dependent && !(init_fun == "zero_init")
        mesh_base = "special_"*init_fun;
    else
        mesh_base = "zero_init";
    end
    
    solvepart = "";
    if prob.time_dependent
        solvepart = "
        //////////////will be generated/////////////////////////////////////////////
        #include \"Stepper.cpp\"
        ////////////////////////////////////////////////////////////////////////////
        ";
    else
        solvepart = "// This uses the generated RHS code to compute the RHS vector.
        rhsVec.computeVec(rhs,rhs,1.0);
        
        // Solve the linear system. 
        lhsMat.cgSolve("*valvec*",rhs,solve_max_iters,solve_tol,0);";
    end
    
    # Just write the whole skeleton
    content = """
    #include "TreeNode.h"
    #include "mpi.h"
    #include "genPts_par.h"
    #include "sfcSort.h"
    #include "mesh.h"
    #include "dendro.h"
    #include "dendroIO.h"
    #include "octUtils.h"
    #include "functional"
    #include "fdCoefficient.h"
    #include "stencil.h"
    #include "rkTransport.h"
    #include "refel.h"
    #include "operators.h"
    #include "cg.h"
    #include <iostream>
    #include <sys/time.h>
    
    #include "linear_skel.h"
    #include "bilinear_skel.h"
    #include "genfunctions.h"

    double rtclock();

    int main (int argc, char** argv){
        
        double t_setup, t_solve;

        MPI_Init(&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;

        int rank, npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &npes);
        
        if(!rank)
            t_setup = rtclock();
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Config.cpp"
        ////////////////////////////////////////////////////////////////////////////
        
        if (argc > 4) {
            m_uiMaxDepth = atoi(argv[1]);
            wavelet_tol = atof(argv[2]);
            partition_tol = atof(argv[3]);
            eOrder = atoi(argv[4]);
        }
        
        Point domain_min(0,0,0);
        Point domain_max(1,1,1);
        
        Point grid_min(0, 0, 0);
        Point grid_max((1u << m_uiMaxDepth), (1u << m_uiMaxDepth), (1u << m_uiMaxDepth));
        
        double Rg_x=(grid_max.x()-grid_min.x());
        double Rg_y=(grid_max.y()-grid_min.y());
        double Rg_z=(grid_max.z()-grid_min.z());

        double Rd_x=(domain_max.x()-domain_min.x());
        double Rd_y=(domain_max.y()-domain_min.y());
        double Rd_z=(domain_max.z()-domain_min.z());

        const Point d_min=domain_min;
        const Point d_max=domain_max;

        const Point g_min=grid_min;
        const Point g_max=grid_max;
        
        std::function<double(double)> gridX_to_X = [d_min,g_min,Rd_x,Rg_x](const double x){
            return d_min.x() + (x-g_min.x())*Rd_x/Rg_x;
        };
        
        std::function<double(double)> gridY_to_Y = [d_min,g_min,Rd_y,Rg_y](const double y){
            return d_min.y() + (y-g_min.y())*Rd_y/Rg_y;
        };
        
        std::function<double(double)> gridZ_to_Z = [d_min,g_min,Rd_z,Rg_z](const double z){
            return d_min.z() + (z-g_min.z())*Rd_z/Rg_z;
        };
        
        std::function<void(double,double,double,double*)> zero_init = [](const double x,const double y,const double z,double *var){
            var[0]=0;
        };
        
        """*functional_defs*"""
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Problem.cpp"
        /////////////////////////////////////////////////////////////////////////////
        
        if (!rank) {
             std::cout << YLW << "maxDepth: " << m_uiMaxDepth << NRM << std::endl;
             std::cout << YLW << "wavelet_tol: " << wavelet_tol << NRM << std::endl;
             std::cout << YLW << "partition_tol: " << partition_tol << NRM << std::endl;
             std::cout << YLW << "eleOrder: " << eOrder << NRM << std::endl;
        }
        
        FinchDendroGenfunctions::setDepth(m_uiMaxDepth);
        
        _InitializeHcurve(m_uiDim);
        RefElement refEl(m_uiDim,eOrder);
        
        ot::DA* octDA=new ot::DA("""* mesh_base *""",1,comm,eOrder,wavelet_tol,100,partition_tol,ot::FEM_CG);
        
        // Variable info will also be generated, but for now assume a single scalar variable
        std::vector<double> uSolVec;
        octDA->createVector(uSolVec,false,false,DOF);
        double *uSolVecPtr=&(*(uSolVec.begin()));

        FinchDendroSkeleton::LHSMat lhsMat(octDA,1);
        lhsMat.setProblemDimensions(domain_min,domain_max);
        lhsMat.setGlobalDofVec(uSolVecPtr);
        
        FinchDendroSkeleton::RHSVec rhsVec(octDA,1);
        rhsVec.setProblemDimensions(domain_min,domain_max);
        rhsVec.setGlobalDofVec(uSolVecPtr);
        
        // This assumes some things
        lhsMat.setBdryFunction(FinchDendroGenfunctions::"""* cpp_gen_string(prob.bc_func[1,1][1]) *""");
        rhsVec.setBdryFunction(FinchDendroGenfunctions::"""* cpp_gen_string(prob.bc_func[1,1][1]) *""");
        
        // Allocate dofs
        """*alloc_dof*"""
        double * rhs=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI_RHS, false,false); // linear part
        
        // Init dofs
        """*set_dof*"""
        octDA->setVectorByFunction(rhs,zero_init,false,false,1); // zeros
        
        // Solve
        if(!rank){
            t_setup = rtclock() - t_setup;
            t_solve = rtclock();
        }
        """*solvepart*"""
        
        if(!rank)
            t_solve = rtclock() - t_solve;
        
        if(!rank){
            std::cout<<"setup time: "<<t_setup<<std::endl;
            std::cout<<"solve time: "<<t_solve<<std::endl;
            std::ofstream logfile;
            logfile.open ("""*"\""*CodeGenerator.genFileName*"_timelog.txt\""*""", std::ios::out | std::ios::app);
            logfile<<"----------------------------------"<<std::endl;
            logfile<<"num. procs: "<<npes<<std::endl;
            logfile<<"setup time: "<<t_setup<<std::endl;
            logfile<<"solve time: "<<t_solve<<std::endl;
            logfile<<"----------------------------------"<<std::endl;
            logfile.close();
        }
        
        // Output
        //////////////will be generated/////////////////////////////////////////////
        #include "Output.cpp"
        ////////////////////////////////////////////////////////////////////////////
        
        octDA->destroyVector(uSolVec);

        if(!rank)
            std::cout<<" End of computation. "<<std::endl;

        delete octDA;

        MPI_Finalize();
        return 0;
    }
    
    double rtclock(){
        struct timezone Tzp;
        struct timeval Tp;
        int stat;
        stat = gettimeofday(&Tp, &Tzp);
        if(stat != 0){
            printf("Error returned from gettimeofday\\n");
        }
        return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
    }
    """
    print(file, content);
end

# This is a class containing various utility functions 
function dendro_utils_file()
    utilsfile = add_generated_file("Utils.cpp", dir="src");
    
    content = "
// Utilities

"
    println(utilsfile, content);
end

function dendro_config_file(dparams=(5, 1, 0.3, 0.000001, 100))
    file = add_generated_file("Config.cpp", dir="src");
    # dparams has (maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
    println(file, "m_uiMaxDepth = "*string(dparams[1])*";           // mesh refinement depth");
    println(file, "double wavelet_tol = "*string(dparams[2])*";     // tolerance for approximating functions(f) determines mesh fineness");
    println(file, "double partition_tol = "*string(dparams[3])*"; // load balancing parameter");
    println(file, "double solve_tol = "*string(dparams[4])*";            // tol used by cgsolve for stopping iterations");
    println(file, "unsigned int solve_max_iters = "*string(dparams[5])*"; // used by cgsolve");
    # From the config struct
    println(file, "unsigned int eOrder  = "*cpp_gen_string(config.basis_order_min)*";");
    
    println(file, "int config_dimension = "*cpp_gen_string(config.dimension)*";"); 
    
    println(file, "const char* config_solver = "*cpp_gen_string(config.solver_type)*";");
    println(file, "const char* config_trial_function = "*cpp_gen_string(config.trial_function)*";");
    println(file, "const char* config_test_function = "*cpp_gen_string(config.test_function)*";");
    println(file, "const char* config_elemental_nodes = "*cpp_gen_string(config.elemental_nodes)*";");
    println(file, "const char* config_quadrature= "*cpp_gen_string(config.quadrature)*";");
end

function dendro_prob_file()
    file = add_generated_file("Problem.cpp", dir="src");
    # variable and coefficient DOFs
    dofs = 0;
    varpart = "";
    varnames = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar)
            varpart = varpart*"M_UI"*string(variables[i].symvar[j]);
            varnames = varnames*"\"m_ui"*string(variables[i].symvar[j])*"\", ";
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
            if prob.time_dependent
                # Also add next versions for the next time step
                varpart = varpart*"M_UI"*string(variables[i].symvar[j])*"_NEXT, ";
                varnames = varnames*"\"m_ui"*string(variables[i].symvar[j])*"_next\", ";
                dofs = dofs + 1;
            end
            dofs = dofs+1;
        end
    end
    coefpart = "";
    coefnames = "";
    # # Now coefficients are not stored as dofs
    # for i=1:length(coefficients)
    #     for j=1:length(coefficients[i].symvar)
    #         coefpart = coefpart*"M_UI"*string(coefficients[i].symvar[j])*", ";
    #         coefnames = coefnames*"\"m_ui"*string(coefficients[i].symvar[j])*"\", ";
    #         dofs = dofs+1;
    #     end
    # end
    dofs += 1; # for linear part
    println(file, "enum VAR{"*varpart*coefpart*"M_UI_RHS}; // variables, coefficients, linear part");
    println(file, "const char * VAR_NAMES[]={"*varnames*coefnames*"\"m_uiRHS\"};");
    println(file, "const unsigned int DOF="*string(dofs)*";"); 
    
    # From the prob struct
    println(file, "const char* prob_bc_type[] = "*cpp_gen_string(prob.bc_type)*";");
    println(file, "int prob_bid[] = "*cpp_gen_string(prob.bid)*";");
    
    # for time dependent problems
    if prob.time_dependent
        println(file, "double tBegin = 0;");
        println(file, "double tEnd = "*string(prob.end_time)*";");
    end
end

# Functions for boundary conditions, coeficients, etc.
function dendro_genfunction_file()
    file = add_generated_file("genfunctions.cpp", dir="src");
    header = add_generated_file("genfunctions.h", dir="include");
    
    content = "
#include \"genfunctions.h\"
#include <math.h>

unsigned int genfunction_maxDepth, genfunction_grid_max;
// set the max depth of the grid
void FinchDendroGenfunctions::setDepth(unsigned int d){
    genfunction_maxDepth = d;
    genfunction_grid_max = 1u<<genfunction_maxDepth;
}
// Converts grid coordinates to spatial coordinates
double FinchDendroGenfunctions::gridCoordToCoord(const double x, const int comp) {
    // cx = cmin + gx * (cmax-cmin) / (gmax-gmin)
    // Note that gmin=0, gmax = grid_max, cmax and cmin are generated.
    double cx=0.0;
    switch(comp){
        case(0):
            cx = 0.0 + x * 1.0 / genfunction_grid_max;
            break;
        case(1):
            cx = 0.0 + x * 1.0 / genfunction_grid_max;
            break;
        case(2):
            cx = 0.0 + x * 1.0 / genfunction_grid_max;
            break;
    }
    return cx;
}

// The genfunctions
";
    
    indent = "";
    args = ["x"; "y"; "z"; "t"];
    argtypes = ["const double"; "const double"; "const double"; "const double"];
    rettype = "double";
    genfunction_list = "";
    for i = 1:length(genfunctions)
        str = dendro_genfunction_to_string(genfunctions[i], xsym=:cx , ysym=:cy , zsym=:cz );
        fun = [];
        try
            num = parse(Float64, str)
            # If this works, just return that number.
            append!(fun, ["return "*str*";"]);
        catch
            append!(fun, [
                "double cx = FinchDendroGenfunctions::gridCoordToCoord(x, 0);",
                "double cy = FinchDendroGenfunctions::gridCoordToCoord(y, 1);",
                "double cz = FinchDendroGenfunctions::gridCoordToCoord(z, 2);",
                "return "*str*";"
                ]);
        end
        
        lines = cpp_function_def(indent, "FinchDendroGenfunctions::"*genfunctions[i].name, args, argtypes, rettype, fun);
        for j=1:length(lines)
            content *= lines[j] * "\n";
        end
        
        genfunction_list *= "    double "*genfunctions[i].name*"(const double x, const double y, const double z, const double t);\n";
    end
    
    print(file, content);
    
    hcontent = "
#ifndef DENDRO_5_0_GENFUNCTIONS_H
#define DENDRO_5_0_GENFUNCTIONS_H

namespace FinchDendroGenfunctions
{
    void setDepth(unsigned int d);
    double gridCoordToCoord(const double x, const int comp);
    
"*genfunction_list*"
}

#endif //DENDRO_5_0_GENFUNCTIONS_H
    ";
    print(header, hcontent);
end

# Time stepper (Only added for time dependent problems)
function dendro_stepper_file()
    file = add_generated_file("Stepper.cpp", dir="src");
    varvec = string(variables[1].symvar[1]);
    varvec_next = string(variables[1].symvar[1])*"_next";
    content = "
double currentT = tBegin;
// Time stepper parameters
int Nsteps = "*string(50)*";
double dt = tEnd / Nsteps;
rhsVec.setDt(dt);
lhsMat.setDt(dt);

for(int ti=0; ti<Nsteps; ti++){
    int activeRank=octDA->getRankActive();
    if(!activeRank)
        std::cout<<\"step \"<<ti<<\" (t=\"<<currentT<<\")\"<<std::endl;
        
    rhsVec.computeVec(rhs,rhs,1.0);
    lhsMat.cgSolve("*varvec_next*",rhs,solve_max_iters,solve_tol,0);
    
    // not sure of the best way to do this. Does this work for now?
    swap("*varvec*", "*varvec_next*");
    
    currentT += dt;
}";
    
    print(file, content);
end

# This is the matrix assembly code.
# It requires code that was generated for this target.
function dendro_bilinear_file(code)
    file = add_generated_file("Bilinear.cpp", dir="src");
    skeleton_file = add_generated_file("bilinear_skel.cpp", dir="src");
    skeleton_headerfile = add_generated_file("bilinear_skel.h", dir="include");
    
    print(file, code);
    
    varpart = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar)
            varpart = varpart*"M_UI"*string(variables[i].symvar[j]);
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
            if prob.time_dependent
                # Also add next versions for the next time step
                varpart = varpart*"M_UI"*string(variables[i].symvar[j])*"_NEXT, ";
            end
        end
    end
    coefpart = "";
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar)
            coefpart = coefpart*"M_UI"*string(coefficients[i].symvar[j])*", ";
        end
    end
    
    content = """
    #include "bilinear_skel.h"

    FinchDendroSkeleton::LHSMat::LHSMat(ot::DA* da,unsigned int dof) : feMatrix(da,dof)
    {
        const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
        imV1=new double[nPe];
        imV2=new double[nPe];

        Qx=new double[nPe];
        Qy=new double[nPe];
        Qz=new double[nPe];
    }

    FinchDendroSkeleton::LHSMat::~LHSMat()
    {
        delete [] imV1;
        delete [] imV2;

        delete [] Qx;
        delete [] Qy;
        delete [] Qz;

        imV1=NULL;
        imV2=NULL;

        Qx=NULL;
        Qy=NULL;
        Qz=NULL;
    }

    void FinchDendroSkeleton::LHSMat::elementalMatVec(const VECType* in,VECType* out, double*coords,double scale)
    {
        const RefElement* refEl=m_uiOctDA->getReferenceElement();

        const double * Q1d=refEl->getQ1d();
        const double * QT1d=refEl->getQT1d();
        const double * Dg=refEl->getDg1d();
        const double * DgT=refEl->getDgT1d();
        const double * W1d=refEl->getWgq();

        const unsigned int eleOrder=refEl->getOrder();
        const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
        const unsigned int nrp=eleOrder+1;

        Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
        Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

        const double refElSz=refEl->getElementSz();
        
        const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
        const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
        const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


        const double Jx = 1.0/(refElSz/(double (szX)));
        const double Jy = 1.0/(refElSz/(double (szY)));
        const double Jz = 1.0/(refElSz/(double (szZ)));
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Bilinear.cpp"
        ////////////////////////////////////////////////////////////////////////////
    }
    
    void FinchDendroSkeleton::LHSMat::setBdryFunction(std::function<double(double,double,double,double)> bdry){
        bdry_function = bdry;
    }
    
    void FinchDendroSkeleton::LHSMat::setGlobalDofVec(double* gdv){
        grandDofVecPtr = gdv;
    }

    bool FinchDendroSkeleton::LHSMat::preMatVec(const VECType* in,VECType* out,double scale)
    {
        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        for(unsigned int i=0;i<bdyIndex.size();i++){
            
            out[bdyIndex[i]] = in[bdyIndex[i]]; // Dirichlet BC
        }

        return true;
    }

    bool FinchDendroSkeleton::LHSMat::postMatVec(const VECType* in,VECType* out,double scale) {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        for(unsigned int i=0;i<bdyIndex.size();i++){
            
            out[bdyIndex[i]] = in[bdyIndex[i]]; // Dirichlet BC
        }

        return true;
    }


    double FinchDendroSkeleton::LHSMat::gridX_to_X(double x)
    {
        double Rg_x=((1u<<m_uiMaxDepth)-0);
        return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
    }

    double FinchDendroSkeleton::LHSMat::gridY_to_Y(double y)
    {
        double Rg_y=((1u<<m_uiMaxDepth)-0);
        return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
    }


    double FinchDendroSkeleton::LHSMat::gridZ_to_Z(double z)
    {
        double Rg_z=((1u<<m_uiMaxDepth)-0);
        return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
    }
    
    void FinchDendroSkeleton::LHSMat::setDt(double newdt)
    {
        dt = newdt;
    }

    int FinchDendroSkeleton::LHSMat::cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var)
    {
        double resid,alpha,beta,rho,rho_1;
        int status=1; // 0 indicates it has solved the system within the specified max_iter, 1 otherwise.

        const unsigned int local_dof=m_uiOctDA->getLocalNodalSz();

        MPI_Comm globalComm=m_uiOctDA->getGlobalComm();

        if(m_uiOctDA->isActive())
        {

            int activeRank=m_uiOctDA->getRankActive();
            int activeNpes=m_uiOctDA->getNpesActive();

            MPI_Comm activeComm=m_uiOctDA->getCommActive();

            double* p;
            double* z;
            double* q;
            double* Ax;
            double* Ap;
            double* r0;
            double* r1;

            m_uiOctDA->createVector(p);
            m_uiOctDA->createVector(z);
            m_uiOctDA->createVector(q);

            m_uiOctDA->createVector(Ax);
            m_uiOctDA->createVector(Ap);
            m_uiOctDA->createVector(r0);
            m_uiOctDA->createVector(r1);

            double normb = normLInfty(b,local_dof,activeComm);
            par::Mpi_Bcast(&normb,1,0,activeComm);

            if(!activeRank)
                std::cout<<"normb = "<<normb<<std::endl;

            matVec(x,Ax);

            /*char fPrefix[256];
            sprintf(fPrefix,"%s_%d","cg",0);
            const char * varNames[]={"U"};
            const double * var[]={Ax};
            io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);
            */
            for(unsigned int i=0;i<local_dof;i++)
            {
                r0[i]=b[i]-Ax[i];
                p[i]=r0[i];
            }

            if (normb == 0.0)
                normb = 1;

            double normr=normLInfty(r0,local_dof,activeComm);
            par::Mpi_Bcast(&normr,1,0,activeComm);
            if(!activeRank) std::cout<<"initial residual : "<<(normr/normb)<<std::endl;

            if ((resid = normr / normb) <= tol) {
                tol = resid;
                max_iter = 0;

                m_uiOctDA->destroyVector(p);
                m_uiOctDA->destroyVector(z);
                m_uiOctDA->destroyVector(q);

                m_uiOctDA->destroyVector(Ax);
                m_uiOctDA->destroyVector(Ap);
                m_uiOctDA->destroyVector(r0);
                m_uiOctDA->destroyVector(r1);

                status=0;
            }

            if(status!=0)
            {
                for(unsigned int i=1;i<=max_iter;i++)
                {
                    matVec(p,Ap);

                    alpha=(dot(r0,r0,local_dof,activeComm)/dot(p,Ap,local_dof,activeComm));
                    par::Mpi_Bcast(&alpha,1,0,activeComm);

                    //if(!activeRank) std::cout<<"rank: " <<activeRank<<" alpha: "<<alpha<<std::endl;
                    for(unsigned int e=0;e<local_dof;e++)
                    {
                        x[e]+=alpha*p[e];
                        r1[e]=r0[e]-alpha*Ap[e];
                    }

                    normr=normLInfty(r1,local_dof,activeComm);
                    par::Mpi_Bcast(&normr,1,0,activeComm);

                    if((!activeRank) && (i%10==0)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;

                    if ((resid = normr / normb) <= tol) {

                        if((!activeRank)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;
                        tol = resid;
                        m_uiOctDA->destroyVector(p);
                        m_uiOctDA->destroyVector(z);
                        m_uiOctDA->destroyVector(q);

                        m_uiOctDA->destroyVector(Ax);
                        m_uiOctDA->destroyVector(Ap);
                        m_uiOctDA->destroyVector(r0);
                        m_uiOctDA->destroyVector(r1);

                        status=0;
                        break;
                    }

                    beta=(dot(r1,r1,local_dof,activeComm)/dot(r0,r0,local_dof,activeComm));
                    par::Mpi_Bcast(&beta,1,0,activeComm);

                    //if(!activeRank) std::cout<<"<r_1,r_1> : "<<dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)<<" <r_0,r_0>: "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" beta "<<beta<<std::endl;

                    for(unsigned int e=0;e<local_dof;e++)
                    {
                        p[e]=r1[e]+beta*p[e];
                        r0[e]=r1[e];
                    }
                }

                if(status!=0)
                {
                    tol = resid;
                    m_uiOctDA->destroyVector(p);
                    m_uiOctDA->destroyVector(z);
                    m_uiOctDA->destroyVector(q);

                    m_uiOctDA->destroyVector(Ax);
                    m_uiOctDA->destroyVector(Ap);
                    m_uiOctDA->destroyVector(r0);
                    m_uiOctDA->destroyVector(r1);
                    status=1;
                }
            }
        }

        // bcast act as a barrier for active and inactive meshes.
        par::Mpi_Bcast(&tol,1,0,globalComm);
        return status;
    }
    """
    print(skeleton_file, content);
    
    content = """
    #ifndef DENDRO_5_0_BILINEAR_SKEL_H
    #define DENDRO_5_0_BILINEAR_SKEL_H

    #include "oda.h"
    #include "feMatrix.h"
    #include "genfunctions.h"

    namespace FinchDendroSkeleton
    {
        class LHSMat : public feMatrix<LHSMat>{

        private:
            // some additional work space variables to perform elemental MatVec
            double* imV1;
            double* imV2;
            double* Qx;
            double* Qy;
            double* Qz;
            // function for boundary
            std::function<double(double,double,double,double)> bdry_function;
            // coefficient vectors
            double* grandDofVecPtr;
            """*"enum VAR{"*varpart*coefpart*"M_UI_RHS};"*"""
            // time step
            double dt;

        public:
            /**@brief: constructor*/
            LHSMat(ot::DA* da,unsigned int dof=1);

            /**@brief default destructor*/
            ~LHSMat();

            /**@biref elemental matvec*/
            virtual void elementalMatVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);
            
            /**@brief set boundary function*/	
            void setBdryFunction(std::function<double(double,double,double,double)> bdry);
            
            /**@brief set pointer to global dof vector*/	
                void setGlobalDofVec(double* gdv);
            
            /**@brief things need to be performed before matvec (i.e. coords transform)*/
            bool preMatVec(const VECType* in,VECType* out,double scale=1.0);

            /**@brief things need to be performed after matvec (i.e. coords transform)*/
            bool postMatVec(const VECType* in,VECType* out,double scale=1.0);

            /**@brief octree grid x to domin x*/
            double gridX_to_X(double x);
            /**@brief octree grid y to domin y*/
            double gridY_to_Y(double y);
            /**@brief octree grid z to domin z*/
            double gridZ_to_Z(double z);
            
            void setDt(double dt);

            int cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var=0);
        };
    }

    #endif //DENDRO_5_0_BILINEAR_SKEL_H
    """
    print(skeleton_headerfile, content);
end

# This is the vector assembly code.
# It requires code that was generated for this target.
function dendro_linear_file(code)
    file = add_generated_file("Linear.cpp", dir="src");
    skeleton_file = add_generated_file("linear_skel.cpp", dir="src");
    skeleton_headerfile = add_generated_file("linear_skel.h", dir="include");
    
    print(file, code);
    
    varpart = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar)
            varpart = varpart*"M_UI"*string(variables[i].symvar[j]);
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
            if prob.time_dependent
                # Also add next versions for the next time step
                varpart = varpart*"M_UI"*string(variables[i].symvar[j])*"_NEXT, ";
            end
        end
    end
    coefpart = "";
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar)
            coefpart = coefpart*"M_UI"*string(coefficients[i].symvar[j])*", ";
        end
    end
    
    content = """
    #include "linear_skel.h"

    FinchDendroSkeleton::RHSVec::RHSVec(ot::DA* da,unsigned int dof) : feVector(da,dof)
    {
        const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
        imV1=new double[nPe];
        imV2=new double[nPe];
    }

    FinchDendroSkeleton::RHSVec::~RHSVec()
    {
        delete [] imV1;
        delete [] imV2;

        imV1=NULL;
        imV2=NULL;
    }

    void FinchDendroSkeleton::RHSVec::elementalComputVec(const VECType* in,VECType* out, double*coords,double scale)
    {
        const RefElement* refEl=m_uiOctDA->getReferenceElement();
        const double * Q1d=refEl->getQ1d();
        const double * QT1d=refEl->getQT1d();
        const double * Dg=refEl->getDg1d();
        const double * W1d=refEl->getWgq();

        const unsigned int eleOrder=refEl->getOrder();
        const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
        const unsigned int nrp=eleOrder+1;

        Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
        Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

        const double refElSz=refEl->getElementSz();
        
        const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
        const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
        const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());
        
        const double Jx = 1.0/(refElSz/(double (szX)));
        const double Jy = 1.0/(refElSz/(double (szY)));
        const double Jz = 1.0/(refElSz/(double (szZ)));
        
        double* Qx = {0};
        double* Qy = {0};
        double* Qz = {0};
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Linear.cpp"
        ////////////////////////////////////////////////////////////////////////////
    }
    
    void FinchDendroSkeleton::RHSVec::setBdryFunction(std::function<double(double,double,double,double)> bdry){
        bdry_function = bdry;
    }
    
    void FinchDendroSkeleton::RHSVec::setGlobalDofVec(double* gdv){
        grandDofVecPtr = gdv;
    }

    bool FinchDendroSkeleton::RHSVec::preComputeVec(const VECType* in,VECType* out, double scale)
    {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);
        
        double x,y,z,t;
        for(unsigned int i=0;i<bdyIndex.size();i++){
            x = bdyCoords[i*3+0];
            y = bdyCoords[i*3+1];
            z = bdyCoords[i*3+2];
            t = 0.0;
            
            out[bdyIndex[i]] = bdry_function(x,y,z,t);
        }

        return true;
    }

    bool FinchDendroSkeleton::RHSVec::postComputeVec(const VECType* in,VECType* out, double scale) {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        double x,y,z,t;
        for(unsigned int i=0;i<bdyIndex.size();i++){
            x = bdyCoords[i*3+0];
            y = bdyCoords[i*3+1];
            z = bdyCoords[i*3+2];
            t = 0.0;
            
            out[bdyIndex[i]] = bdry_function(x,y,z,t);
        }

        return true;
    }


    double FinchDendroSkeleton::RHSVec::gridX_to_X(double x)
    {
        double Rg_x=((1u<<m_uiMaxDepth)-0);
        return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
    }

    double FinchDendroSkeleton::RHSVec::gridY_to_Y(double y)
    {
        double Rg_y=((1u<<m_uiMaxDepth)-0);
        return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
    }


    double FinchDendroSkeleton::RHSVec::gridZ_to_Z(double z)
    {
        double Rg_z=((1u<<m_uiMaxDepth)-0);
        return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
    }
    
    void FinchDendroSkeleton::RHSVec::setDt(double newdt)
    {
        dt = newdt;
    }
    
    void FinchDendroSkeleton::RHSVec::getCorrectElementalValues(const VECType *in, VECType* out) {
        VECType* _in=NULL;
        
        m_uiOctDA->nodalVecToGhostedNodal(in,_in,false,m_uiDof);
        
        m_uiOctDA->readFromGhostBegin(_in, m_uiDof);
        
        m_uiOctDA->getElementNodalValues(_in, out, m_uiOctDA->curr(), m_uiDof);
        
        m_uiOctDA->readFromGhostEnd(_in, m_uiDof);
        
        m_uiOctDA->destroyVector(_in);
        
        return;
    }
    """
    print(skeleton_file, content);
    
    content = """
    #ifndef DENDRO_5_0_LINEAR_SKEL_H
    #define DENDRO_5_0_LINEAR_SKEL_H

    #include "oda.h"
    #include "feVector.h"
    #include "genfunctions.h"

    namespace FinchDendroSkeleton
    {
        class RHSVec : public feVector<RHSVec>{

        private:

            double * imV1;
            double * imV2;
            // function for boundary
            std::function<double(double,double,double,double)> bdry_function;
            // coefficient vectors
            double* grandDofVecPtr;
            """*"enum VAR{"*varpart*coefpart*"M_UI_RHS};"*"""
            // time step
            double dt;

        public:
            RHSVec(ot::DA* da,unsigned int dof=1);
            ~RHSVec();
            
            /**@biref elemental compute vec for rhs*/
            virtual void elementalComputVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);
            
            /**@brief set boundary function*/	
            void setBdryFunction(std::function<double(double,double,double,double)> bdry);
            
            /**@brief set pointer to global dof vector*/	
                void setGlobalDofVec(double* gdv);
            
            bool preComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            bool postComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            /**@brief octree grid x to domin x*/
            double gridX_to_X(double x);
            /**@brief octree grid y to domin y*/
            double gridY_to_Y(double y);
            /**@brief octree grid z to domin z*/
            double gridZ_to_Z(double z);
            
            void setDt(double dt);
            
            void getCorrectElementalValues(const VECType *in, VECType* out);
        };
    }

    #endif //DENDRO_5_0_LINEAR_SKEL_H
    """
    print(skeleton_headerfile, content);
end

# Right now this just tries to plot
function dendro_output_file()
    file = add_generated_file("Output.cpp", dir="src");
    
    if config.output_format == VTK
        println(file, "octDA->vecTopvtu(uSolVecPtr,\""*CodeGenerator.genFileName*"\",(char**)VAR_NAMES,false,false,DOF);");
    end
end

function dendro_cmake_file()
    file = add_generated_file("CMakeLists.txt", dir="", make_header_text=false);
    content = """
cmake_minimum_required(VERSION 2.8)
project("""*CodeGenerator.genFileName*""")

set("""*CodeGenerator.genFileName*"""_INC include/linear_skel.h
            include/bilinear_skel.h
            include/genfunctions.h
        )

set("""*CodeGenerator.genFileName*"""_SRC src/linear_skel.cpp
            src/bilinear_skel.cpp
            src/genfunctions.cpp
        )

set(SOURCE_FILES src/"""*CodeGenerator.genFileName*""".cpp \${"""*CodeGenerator.genFileName*"""_INC} \${"""*CodeGenerator.genFileName*"""_SRC})
add_executable("""*CodeGenerator.genFileName*""" \${SOURCE_FILES})
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/include/test)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${MPI_INCLUDE_PATH})
target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories("""*CodeGenerator.genFileName*""" PRIVATE \${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries("""*CodeGenerator.genFileName*""" dendro5 \${LAPACK_LIBRARIES} \${MPI_LIBRARIES} m)

""";
    print(file, content);
end

function dendro_readme_file()
    file = add_generated_file("README.txt", dir="", make_header_text=false);
    content = """
Basic instructions for compiling this generated code with dendro.
    1. Place this generated directory in the Dendro directory.
    
    2. Append the following line to the Dendro CMakeLists.txt file:
        add_subdirectory("""*uppercasefirst(CodeGenerator.genFileName)*""")
    
    3. Remake Dendro as usual. For example, to build to the directory "build",
       enter the Dendro directory and do the following.
        a. If it doesn't exist, make the build directory: "mkdir build"
        b. "cd build"
        c. "ccmake .."
        d. "make"
        e. "cd """*uppercasefirst(CodeGenerator.genFileName)*""" "
    
    4. Run the executable with "./"""*CodeGenerator.genFileName*""" "
    
    *Note: The directory will have an upper case first letter. 
           The executable will match the name passed to generateFor().
"""
    print(file, content);
end

###########################################################################################################
# Symbolic to code layer generation

  ######      #####     ######    #######       ##           ###     ##    ##  #######  ######
##     ##   ###   ###   ##   ##   ##            ##          ## ##     ##  ##   ##       ##   ##
##          ##     ##   ##    ##  ######        ##         ##   ##     ####    ######   ######
##     ##   ###   ###   ##   ##   ##            ##        #########     ##     ##       ##  ##
  ######      #####     ######    #######       #######  ##       ##    ##     #######  ##   ##

###########################################################################################################

# Allocate, compute, or fetch all needed values
function dendrotarget_prepare_needed_values(entities, var, lorr, vors)
    to_delete = []; # arrays allocated with new that need deletion
    
    code = "";
    coef_loop = "";
    coef_interp = "";
    for i=1:length(entities)
        cname = CodeGenerator.make_coef_name(entities[i]);
        if CodeGenerator.is_test_function(entities[i])
            # Assign it a transpose quadrature matrix
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= "// " * cname * " = TRQ"*string(entities[i].derivs[di])*"; // d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                end
            else
                code *= "// " * cname * " = refel.Q'; // test function.\n";
            end
        elseif CodeGenerator.is_unknown_var(entities[i], var) && lorr == LHS
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= "// " * cname * " = RQ"*string(entities[i].derivs[di])*"; // d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                end
            else
                code *= "// " * cname * " = refel.Q; // trial function.\n";
            end
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = CodeGenerator.get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt
            elseif ctype == 0
                # It was a number, do nothing?
            elseif ctype == 1 # a constant wrapped in a coefficient
                # This generates something like: coef_k_i = 4;
                if length(entities[i].derivs) > 0
                    code *= "double " * cname * " = 0.0; // NOTE: derivative applied to constant coefficient = 0\n";
                else
                    code *= "double " * cname * " = " * string(cval) * ";\n";
                end
                
            elseif ctype == 2 # a coefficient function
                if vors == "volume"
                    push!(to_delete, cname);
                    short_name = split(cname, "coef_")[end];
                    code *= "double* "*cname*" = new double[nPe];\n";
                    # code *= "m_uiOctDA->getElementNodalValues(m_uiOctDA->getVecPointerToDof(grandDofVecPtr, VAR::M_UI"*short_name*", false,false), "*cname*", m_uiOctDA->curr(), m_uiDof);\n";
                    # BREAKS?? -> code *= "FinchDendroSkeleton::RHSVec::getCorrectElementalValues(m_uiOctDA->getVecPointerToDof(grandDofVecPtr, VAR::M_UI"*short_name*", false,false), "*cname*");\n";
                    coef_loop *= "    "*cname*"[coefi] = FinchDendroGenfunctions::genfunction_"*string(cval-1)*"(coords[coefi*3+0], coords[coefi*3+1], coords[coefi*3+2], 0.0);\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        println("coefficient derivatives not ready: "*string(entities[i]))
                        # xyzchar = ["x","y","z"];
                        # for di=1:length(entities[i].derivs)
                        #     code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                        #             "; // Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        # end
                    else
                        if lorr==LHS
                            coef_interp *=
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,"*cname*",imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*cname*");\n\n";
                        end
                    end
                else
                    #TODO surface
                end
                
            elseif ctype == 3 # a known variable value
                # This should only occur for time dependent problems
                # Use the values from the previous time step
                if vors == "volume"
                    push!(to_delete, cname);
                    short_name = split(cname, "coef_")[end];
                    code *= "double* "*cname*" = new double[nPe];\n";
                    code *= "m_uiOctDA->getElementNodalValues(m_uiOctDA->getVecPointerToDof(grandDofVecPtr, VAR::M_UI"*short_name*", false,false), "*cname*", m_uiOctDA->curr(), m_uiDof);\n";
                    # code *= "FinchDendroSkeleton::RHSVec::getCorrectElementalValues(m_uiOctDA->getVecPointerToDof(grandDofVecPtr, VAR::M_UI"*short_name*", false,false), "*cname*");\n";
                else
                    #TODO surface
                end
                
            end
        end # if coefficient
    end # entity loop
    
    # Compute any coefficients at nodes
    code *= 
"for(int coefi=0; coefi<nPe; coefi++){
"*coef_loop*"
}";
    
    if length(coef_interp)>0
        code *= "// Interpolate at quadrature points.\n";
        code *= coef_interp;
    end
    
    return (code, to_delete);
end

function dendrotarget_make_elemental_computation(terms, to_delete, var, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    
    # find some useful numbers
    dofsper = 0;
    offset_ind = [0];
    if typeof(var) <:Array
        offset_ind = zeros(Int, varcount);
        dofsper = length(var[1].symvar);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + length(var[i].symvar);
        end
    else
        dofsper = length(var.symvar);
    end
    
    # This will hold the full code string
    code = "";
    # If temp storage was needed for multiple terms, allocate
    temp_alloc = "";
    dealloc = "";
    # The pieces of the code
    preloop_code = "";
    loop_code = "";
    postloop_code = "";
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        printerr("dendro not ready for multi dofs per node. ");
        # Submatrices or subvectors for each component
        if lorr == LHS
            submatrices = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            submatrices = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(submatrices)
            submatrices[smi] = "";
        end
        
        if typeof(var) <: Array
            term_number = 0;
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        term_number += 1;
                        (trial_code, coef_code, test_code, test_ind, trial_ind) = dendrotarget_generate_term_calculation(terms[vi][ci][i], term_number, var, lorr);
                        temp_alloc *= "double* out_"*string(term_number)*" = new double[nPe];\n";
                        push!(to_delete, "out_"*string(term_number));
                        
                        # println(terms)
                        # println(terms[vi])
                        # println(terms[vi][ci])
                        # println(terms[vi][ci][i])
                        # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
                        # Find the appropriate submatrix for this term
                        submati = offset_ind[vi] + test_ind;
                        submatj = trial_ind;
                        if lorr == LHS
                            submat_ind = submati + dofsper * (submatj-1);
                        else
                            submat_ind = submati;
                        end
                        
                        preloop_code *= trial_code * "\n";
                        loop_code *= coef_code * "\n";
                        postloop_code *= test_code * "\n";
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    term_number += 1;
                    (trial_code, coef_code, test_code, test_ind, trial_ind) = dendrotarget_generate_term_calculation(terms[ci][i], term_number, var, lorr);
                    temp_alloc *= "double* out_"*string(term_number)*" = new double[nPe];\n";
                    push!(to_delete, "out_"*string(term_number));
                    
                    # Find the appropriate submatrix for this term
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                    else
                        submat_ind = test_ind;
                    end
                    
                    preloop_code *= trial_code * "\n";
                    loop_code *= coef_code * "\n";
                    postloop_code *= test_code * "\n";
                end
            end
            
        end
        
    else # one dof
        terms = terms[1];
        #process each term
        for i=1:length(terms)
            (trial_code, coef_code, test_code, test_ind, trial_ind) = dendrotarget_generate_term_calculation(terms[i], i, var, lorr);
            temp_alloc *= "double* out_"*string(i)*" = new double[nPe];\n";
            push!(to_delete, "out_"*string(i));
            
            preloop_code *= string(trial_code) * "\n";
            loop_code *= string(coef_code) * "\n";
            postloop_code *= string(test_code) * "\n";
        end
    end
    
    for i=1:length(to_delete)
        dealloc *= "delete [] "*to_delete[i]*";\n";
    end
    
    # Put the pieces together
    code = "// Allocate temporary storage for each term.\n"
    code *= temp_alloc;
    if lorr == RHS && length(loop_code) > 0
        code *= "double* rhscoefvec = new double[nPe];\n";
        dealloc *= "delete [] rhscoefvec;\n";
    end
    code *= "\n";
    
    if lorr==LHS
        code *= "// Trial function factors.\n"
    end
    code *= preloop_code;
    code *= "// Multiply quadrature weights, geometric factors and coefficients.\n"
    code *= 
"unsigned int node_index;
for(unsigned int k=0;k<(eleOrder+1);k++){
    for(unsigned int j=0;j<(eleOrder+1);j++){
        for(unsigned int i=0;i<(eleOrder+1);i++){
            node_index = (k*nrp+j)*nrp+i;
";
    code *= loop_code;
    code *= "
        }
    }
}
\n";
    code *= "// multiply by test functions\n"
    code *= postloop_code;
    
    combine_loop = "out_1[i]";
    for i=2:length(terms)
        combine_loop *= "+out_"*string(i)*"[i]";
    end
    code *= 
"
// Combine all of the terms into the output.
for(unsigned int i=0;i<nPe;i++){
    out[i]="*combine_loop*";
}\n";
    
    code *= dealloc;
    
    return code;
end

function dendrotarget_generate_term_calculation(term, thisterm, var, lorr)
    test_code = "";
    trial_code = "";
    coef_code = "";
    test_deriv = 0;
    trial_deriv = 0;
    
    # The output of this term will be stored in a temporary array "out_n"
    out_name = "out_"*string(thisterm);
    
    # extract each of the factors.
    (test_part, trial_part, coef_part, test_ind, trial_ind) = CodeGenerator.separate_factors(term, var);
    
    # strip off all negatives, combine and apply the result somewhere convenient
    neg = false;
    if typeof(test_part) == Expr && test_part.args[1] === :- && length(test_part.args) == 2
        neg = !neg;
        test_part = test_part.args[2];
    end
    if typeof(trial_part) == Expr && trial_part.args[1] === :- && length(trial_part.args) == 2
        neg = !neg;
        trial_part = trial_part.args[2];
    end
    if typeof(coef_part) == Expr && coef_part.args[1] === :- && length(coef_part.args) == 2
        neg = !neg;
        coef_part = coef_part.args[2];
    end
    
    # Turn each part into a code string
    # First the test function
    if !(test_part === nothing)
        test_code = "// Multiply by test function factor: "*string(test_part)*"\n";
        if length(test_part.derivs) > 0
            test_deriv = test_part.derivs[1];
            if test_part.derivs[1] == 1
                test_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,"*out_name*",imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,"*out_name*");\n";
            elseif test_part.derivs[1] == 2
                test_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,"*out_name*",imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,"*out_name*");\n";
            elseif test_part.derivs[1] == 3
                test_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,"*out_name*",imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,"*out_name*");\n";
            elseif test_part.derivs[1] == 4
                # This will eventually be ??
                printerr("Derivative index problem in "*string(test_part));
            else
                printerr("Derivative index problem in "*string(test_part));
            end
        else # no derivatives
            test_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,"*out_name*",imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,"*out_name*");\n";
        end
    else
        # There was no test function on this term. This shouldn't happen
        printerr("No test function found for term: "*string(term)*", expect an error.");
        return (-1, -1, -1, -1, -1, -1);
    end
    
    # Then the trial function
    if !(trial_part === nothing) && lorr == LHS
        trial_code = "// Multiply by trial function factor: "*string(trial_part)*"\n";
        if length(trial_part.derivs) > 0
            trial_deriv = trial_part.derivs[1];
            if trial_part.derivs[1] == 1
                trial_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
            elseif trial_part.derivs[1] == 2
                trial_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
            elseif trial_part.derivs[1] == 3
                trial_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,"*out_name*");\n";
            elseif trial_part.derivs[1] == 4
                # This will eventually be ??
                printerr("Derivative index problem in "*string(trial_part));
            else
                printerr("Derivative index problem in "*string(trial_part));
            end
        else # no derivatives
            trial_code *= 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
        end
    else
        # There was no trial function on this term.
    end
    
    # build weight/coefficient parts
    idxstr = "[i]";
    weight_code = "W1d[i]";
    if config.dimension == 2
        idxstr = "[j*nrp+i]";
        weight_code = "W1d[i]*W1d[j]";
    elseif config.dimension == 3
        idxstr = "[(k*nrp+j)*nrp+i]";
        weight_code = "W1d[i]*W1d[j]*W1d[k]";
    end
    
    # If term is negative, apply it here
    if neg
        weight_code = "-"*weight_code;
    end
    
    # Multiply by coefficients inside integral if needed
    # need to change the index for rhs
    if lorr == RHS
        idxstr = "coefi";
    else
        idxstr = "node_index";
    end
    if !(coef_part === nothing) || (lorr == RHS && !(trial_part === nothing))
        if !(coef_part === nothing)
            coef_string = string(dendro_change_math_ops(CodeGenerator.replace_entities_with_symbols(coef_part, index=idxstr)));
        else
            coef_string = "";
        end
        # for time dependent parts there can be variable parts on RHS
        if lorr == RHS && !(trial_part === nothing)
            if length(coef_string)>0
                coef_string = coef_string * " * " * string(dendro_change_math_ops(CodeGenerator.replace_entities_with_symbols(trial_part, index=idxstr)));
            else
                coef_string = string(dendro_change_math_ops(CodeGenerator.replace_entities_with_symbols(trial_part, index=idxstr)));
            end
        end
        
        if lorr == LHS
            weight_code = weight_code*" * "*coef_string;
        else
            rhscoef_code = 
"for(int coefi=0; coefi<nPe; coefi++){
    rhscoefvec[coefi] = "*coef_string*";
}
";
        end
        
    end
    
    # build the inner weight/coef loop
    idxstr = "[node_index]";
    jacobian_factors = ["Jx*Jy*Jz"  "Jy*Jz"     "Jx*Jz"     "Jx*Jy";
                        "Jy*Jz"     "Jy*Jz/Jx"  "Jz"        "Jy";
                        "Jx*Jz"     "Jz"        "Jx*Jz/Jy"  "Jy";
                        "Jx*Jy"     "Jy"        "Jx"        "Jx*Jy/Jz"];
    coef_code = out_name*idxstr*"*=(("*jacobian_factors[test_deriv+1, trial_deriv+1]*")*"*weight_code*");";
    
    # If there was no trial part, it's an RHS and we need to finish the quadrature with this
    if lorr == RHS
        trial_code = rhscoef_code*
"// For RHS, multiply by refel.Q to interpolate to quadrature points.
DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,rhscoefvec,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
    end
    
    return (trial_code, coef_code, test_code, test_ind, trial_ind);
end
