#=
This contains all of the pieces needed to add a new code gen target.
The following three functions must be provided.

1. get_external_language_elements() - file extensions, comment chars etc.
2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code

You will have access to anything in the Finch module scope because this file will be included there.
=#

function get_external_language_elements()
    file_extension = ".m";
    comment_char = "%";
    block_comment = ["%{"; "%}"];
    
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
    
    # Determine if derivative matrices will be required
    need_derivs = false;
    for i=1:length(entities)
        if length(entities[i].derivs) > 0
            need_derivs = true;
            break;
        end
    end
    if need_derivs
        code *= matlabtarget_build_derivative_matrices(lorr, vors);
        code *= "\n";
    end
    
    # Allocate/Evaluate or fetch the values for each needed entity.
    # Note: This could be staged and interleaved with the calculation part for complex problems. TODO
    code *= matlabtarget_prepare_needed_values(entities, var, lorr, vors);
    code *= "\n";
    
    # Form the final elemental calculation
    code *= matlabtarget_make_elemental_computation(terms, var, lorr, vors);
    
    return code;
end

# Create all of the code files.
# The input are code strings generated above. They will be inserted into their appropriate files.
function generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf)
    # Again, this can be split up into smaller functions as you wish.
    matlab_utils_file();
    matlab_main_file(var);
    matlab_config_file();
    matlab_prob_file(var);
    matlab_mesh_file();
    matlab_genfunction_file();
    matlab_bilinear_file(lhs_vol, var);
    matlab_linear_file(rhs_vol, var);
    # matlab_stepper_file(); # TODO
    matlab_output_file();
end

#########################################################
# code writing utilities
#########################################################

# numbers: 2 -> "2"
# strings: "thing" -> "'thing'"
# arrays: [1 2; 3 4] -> "[1 2; 3 4]"
function matlab_gen_string(v)
    if typeof(v) == String
        return "'"*v*"'";
        
    elseif typeof(v) <: Number
        return string(v);
        
    elseif typeof(v) <: Array
        if length(v)>0 && typeof(v[1]) == String
            # need to use a cell array
            bracket1 = "{";
            bracket2 = "}";
        else
            bracket1 = "[";
            bracket2 = "]";
        end
        if ndims(v) == 1 # "[a; b; c]"
            n = length(v);
            str = bracket1;
            for i=1:n
                str = str*matlab_gen_string(v[i]);
                if i < n
                    str = str*"; ";
                end
            end
            str = str*bracket2;
        elseif ndims(v) == 2 # "[a b ; c d]
            (n,m) = size(v);
            str = bracket1;
            for i=1:n
                for j=1:m
                    str = str*matlab_gen_string(v[i,j])*" ";
                end
                if i < n
                    str = str*"; ";
                end
            end
            str = str*bracket2;
        end
        return str;
        
    elseif typeof(v) == GenFunction
        return v.name;
    else
        return string(v);
    end
end

# Returns a string like "fread(f, [3, 15], 'double')" for an array with size (3,15) and type Float64
function matlab_fread(A)
    if typeof(A) <: Array
        if length(A) == 0
            return "0";
        end
        if typeof(A[1]) == Int
            typ = "\'int64\'";
        else
            typ = "\'double\'";
        end
        if length(size(A)) == 1
            sz = "["*string(length(A))*",1]";
        elseif length(size(A)) == 2
            sz = "["*string(size(A,1))*","*string(size(A,2))*"]";
        else # Matlab's fread can only handle 2D arrays! Just smash the rest together.
            sz = "[" * string(size(A,1));
            second_dim = 1;
            for i=2:length(size(A))
                second_dim = second_dim * size(A,i);
            end
            sz *= ", " * string(second_dim) * "]";
        end
        return "fread(f, "*sz*", "*typ*")"
    else
        if typeof(A) == Int64
            typ = "\'int64\'";
        else
            typ = "\'double\'";
        end
        return "fread(f, [1], "*typ*")"
    end
end

# produces code to read a binary struct into matlab
# The struct is labeled with the name and has the same fieldnames as s
function matlab_struct_reader(name, s)
    code = "";
    for fn in fieldnames(typeof(s))
        f = getfield(s,fn);
        if typeof(f) <: Array && length(f) > 0 && typeof(f[1]) <: Array
            code = code * name * "." * string(fn) * " = cell([" * string(length(f)) * ", 1]);\n";
            for i=1:length(f)
                code = code * name * "." * string(fn) * "{" * string(i) * "} = " * matlab_fread(f[i]) * ";\n";
            end
        else
            code = code * name * "." * string(fn) * " = " * matlab_fread(f) * ";\n";
        end
    end
    return code;
end

#######################################################
# Write code files
#######################################################

function matlab_main_file(var)
    file = add_generated_file(project_name*".m", dir="src");
    
    println(file, "clear;");
    println(file, "");
    # These are always included
    println(file, "");
    println(file, "Utils;");
    println(file, "Config;");
    println(file, "Mesh;");
    println(file, "Genfunction;");
    println(file, "Problem;");
    println(file, "Bilinear;");
    println(file, "Linear;");
    println(file, "");
    
    # No time stepping yet
    # Just a simple \ linear solve
    println(file, "result = LHS\\RHS;");
    
    # place values in named variable arrays
    if typeof(var) <: Array
        counter = 0;
        for vi=1:length(var)
            for i=1:var[vi].total_components
                counter += 1;
                println(file, string(var.symbol) * "_"*string(counter)*" = result("*string(counter)*":dofspernode:length(result));");
            end
        end
        println(file, "total_components = "*string(counter)*";");
    elseif var.total_components > 1
        for i=1:var.total_components
            println(file, string(var.symbol) * "_"*string(i)*" = result("*string(i)*":dofspernode:length(result));");
        end
        println(file, "total_components = "*string(var.total_components)*";");
    else
        # single component
        println(file, string(var.symbol) * " = result;");
        println(file, "total_components = 1;");
    end
    
    println(file, "");
    
    println(file, "Output;");
    
end

# This is a class containing various utility functions like tensor math, geometric factors, linear solve
function matlab_utils_file()
    utilsfile = add_generated_file("Utils.m", dir="src");
    
    content = "
classdef Utils
    
    methods(Static)
        % 2D routines
        function y = tensor_IAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N);
            y = y(:);
        end

        function y = tensor_AIX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N)';
            y = y'; 
            y = y(:);
        end

        % 3D routines
        function y = tensor_IIAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N*N);
            y = y(:);
        end

        function y = tensor_IAIX (A, x)
            N = size (A, 1);
            q = reshape(x, N, N, N);
            y = zeros(N,N,N);
            for i=1:N
                y(i,:,:) = A * squeeze( q(i,:,:) );
            end
            y = y(:);
        end

        function y = tensor_AIIX (A, x)
            N = size (A, 1);
            y = reshape(x, N*N, N) * A';
            y = y(:);
        end

        function du = tensor_grad(refel, u)
            du = zeros(length(u), refel.dim);
            if (refel.dim == 2)
                du(:,1) = Utils.tensor_IAX (refel.Dr, u);
                du(:,2) = Utils.tensor_AIX (refel.Dr, u);
            else
                du(:,1) = Utils.tensor_IIAX (refel.Dr, u);
                du(:,2) = Utils.tensor_IAIX (refel.Dr, u);
                du(:,3) = Utils.tensor_AIIX (refel.Dr, u);
            end
        end

        function [dx, dy] = tensor_grad2(A, x)
            dx = Utils.tensor_IAX (A, x);
            dy = Utils.tensor_AIX (A, x);
        end

        function [dx, dy, dz] = tensor_grad3(A, x)
            dx = Utils.tensor_IIAX (A, x);
            dy = Utils.tensor_IAIX (A, x);
            dz = Utils.tensor_AIIX (A, x);
        end


        function [J, D] = geometric_factors(refel, pts)

            if (refel.dim == 0)
                xr  = [1];
                J = xr;
            elseif (refel.dim == 1)
                xr  = refel.Dr*pts;
                J = xr;
            elseif (refel.dim == 2)
                if refel.Nfaces == 3 % triangle
                    xr = refel.Ddr*pts(1,:)';
                    xs = refel.Dds*pts(1,:)';
                    yr = refel.Ddr*pts(2,:)';
                    ys = refel.Dds*pts(2,:)';
                    J = -xs.*yr + xr.*ys;
                    J = J(1); 
                else % quad
                    [xr, xs] = Utils.tensor_grad2 (refel.Dg, pts(1,:));
                    [yr, ys] = Utils.tensor_grad2 (refel.Dg, pts(2,:));

                    J = -xs.*yr + xr.*ys;
                end

            else
                [xr, xs, xt] = Utils.tensor_grad3 (refel.Dg, pts(1,:));
                [yr, ys, yt] = Utils.tensor_grad3 (refel.Dg, pts(2,:));
                [zr, zs, zt] = Utils.tensor_grad3 (refel.Dg, pts(3,:));

                J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
            end

            if (nargout > 1)
                if (refel.dim == 1)
                    D.rx = 1./J;
                elseif (refel.dim == 2)
                    if refel.Nfaces == 3 % triangle
                        D.rx =  ys./J;
                        D.sx = -yr./J;
                        D.ry = -xs./J;
                        D.sy =  xr./J;
                    else % quad
                        D.rx =  ys./J;
                        D.sx = -yr./J;
                        D.ry = -xs./J;
                        D.sy =  xr./J;
                    end

                else
                    D.rx =  (ys.*zt - zs.*yt)./J;
                    D.ry = -(xs.*zt - zs.*xt)./J;
                    D.rz =  (xs.*yt - ys.*xt)./J;

                    D.sx = -(yr.*zt - zr.*yt)./J;
                    D.sy =  (xr.*zt - zr.*xt)./J;
                    D.sz = -(xr.*yt - yr.*xt)./J;

                    D.tx =  (yr.*zs - zr.*ys)./J;
                    D.ty = -(xr.*zs - zr.*xs)./J;
                    D.tz =  (xr.*ys - yr.*xs)./J;
                end

            end
        end
        
        % pack and unpack variables into a global vector
        % packed vars are in a vector, unpacked are in a cell table
        function pk = pack_vars(upk, pk)
            for i=1:Nvars
                pk(i:Nvars:length(pk)) = upk{i};
            end
        end
        
        function upk = unpack_vars(pk)
            upk = cell(Nvars,1);
            for i=1:Nvars
                upk{i} = pk(i:Nvars:length(pk));
            end
        end
    end
end

"
    println(utilsfile, content);
end

function matlab_config_file()
    file = add_generated_file("Config.m", dir="src");
    # Duplicate the config struct
    for f in fieldnames(Finch.Finch_config)
        println(file, "config."*string(f)*" = "*matlab_gen_string(getfield(Finch.config, f))*";");
    end
    println(file, "order = config.basis_order_min;");
end

function matlab_prob_file(var)
    file = add_generated_file("Problem.m", dir="src");
    # Duplicate the prob struct
    for f in fieldnames(Finch.Finch_prob)
        println(file, "prob."*string(f)*" = "*matlab_gen_string(getfield(Finch.prob, f))*";");
    end
    # count dofs per node
    multivar = typeof(var) <:Array;
    dofsper = 0;
    if multivar
        dofsper = var[1].total_components;
        for i=2:length(var)
            dofsper = dofsper + var[i].total_components;
        end
    else
        dofsper = var.total_components;
    end
    println(file, "dofspernode = "*string(dofsper)*";");
end

#=
The mesh file contains any code related to setting up the mesh.
The meshdata file contains all of the data from the Refel, MeshData and Grid structs
These are to be read into matlab by a matlab matlab function in the mesh file.
=#
function matlab_mesh_file()
    file = add_generated_file("Mesh.m", dir="src");
    
    println(file, "f = fopen('MeshData','r');");
    println(file, "% Reference element");
    println(file, matlab_struct_reader("refel", refel));
    println(file, "% mesh data");
    println(file, matlab_struct_reader("mesh_data", mesh_data));
    println(file, "% grid data");
    println(file, matlab_struct_reader("grid_data", grid_data));
    println(file, "fclose(f);");
    
    datfile = add_generated_file("MeshData", dir="src", make_header_text=false);
    # refel
    Finch.CodeGenerator.write_binary_struct(datfile, refel);
    # mesh_data
    Finch.CodeGenerator.write_binary_struct(datfile, mesh_data);
    # grid_data
    Finch.CodeGenerator.write_binary_struct(datfile, grid_data);
end

# Functions for boundary conditions, coeficients, etc.
function matlab_genfunction_file()
    file = add_generated_file("Genfunction.m", dir="src");
    
    for i = 1:length(genfunctions)
        println(file, genfunctions[i].name*"_fun = @(x,y,z,t) ("*genfunctions[i].str*");");
        # Evaluate them at grid points and make them vectors. Make sense???
        println(file, genfunctions[i].name*" = evaluate_genfun("*genfunctions[i].name*"_fun, grid_data.allnodes, 0);");
    end
    
    # assign variable and coefficient symbols to these vectors
    nvars = 0
    for v in variables
        println(file, string(v.symbol)*" = '"*string(v.symbol)*"';");
        nvars += size(v.values,1);
    end
    println(file, "Nvars = " * string(nvars) * ";");
    for v in coefficients
        if typeof(v.value[1]) == GenFunction
            println(file, string(v.symbol)*" = "*v.value[1].name*";");
        else
            println(file, string(v.symbol)*" = "*string(v.value[1])*";");
        end
    end
    
    evalgenfun = 
"""
function u = evaluate_genfun(genfun, pts, t)
    n = size(pts,2);
    dim = size(pts,1);
    u = zeros(n,1);
    x = 0;
    y = 0;
    z = 0;
    for i=1:n
        x = pts(1,i);
        if dim > 1
            y = pts(2,i);
            if dim > 2
                z = pts(3,i);
            end
        end
        u(i) = genfun(x,y,z,t);
    end
end
"""
    println(file, evalgenfun);
end

# This is the matrix assembly loop.
# It requires code that was generated for this target.
function matlab_bilinear_file(code, var)
    file = add_generated_file("Bilinear.m", dir="src");
    
    Ddotn_part = "
                                Ddotn = RD1(bdrylocal(j),:) * normal(1)";
                                
    if config.dimension > 1
        Ddotn_part *= " + ...
                                RD2(bdrylocal(j),:) * normal(2)";
    end
    if config.dimension > 2
        Ddotn_part *= " + ...
                                RD3(bdrylocal(j),:) * normal(3)";
    end
    Ddotn_part *= ";";
    
    # insert the code part into this skeleton
    content = 
"
dof = size(grid_data.allnodes,2);
ne  = mesh_data.nel;
Np = refel.Np;
glbdof = zeros(Np*dofspernode,1);
I = zeros(ne * Np*Np*dofspernode*dofspernode, 1);
J = zeros(ne * Np*Np*dofspernode*dofspernode, 1);
val = zeros(ne * Np*Np*dofspernode*dofspernode, 1);

%temporary
time = 0;
dt = 1;
var = 0;

% loop over elements
for e=1:ne
    glb = grid_data.loc2glb(:,e)\';
    nodex = grid_data.allnodes(:,glb);
    [detJ, Jac]  = Utils.geometric_factors(refel, nodex);
    wdetj = refel.wg .* detJ;
    
    for di=1:dofspernode
        glbdof(((di-1)*Np + 1):(di*Np)) = (glb - 1).*dofspernode + di;
    end
    
    ind1 = repmat(glbdof,Np*dofspernode,1);
    ind2 = reshape(repmat(glbdof',Np*dofspernode,1),Np*Np*dofspernode*dofspernode,1);
    
    st = (e-1)*Np*Np*dofspernode*dofspernode+1;
    en = e*Np*Np*dofspernode*dofspernode;
    I(st:en) = ind1;
    J(st:en) = ind2;
    
    % Generated elemental computation ############################
"*code*"
    % End generated elemental computation ########################
    
    % handle boundary conditions
    for f=1:length(grid_data.element2face(:,e))
        fid = grid_data.element2face(f,e);
        bid = grid_data.facebid(fid);
        if bid > 0
            % boundary face
            for di=1:dofspernode
                bdrylocaldof = refel.face2local{f} + Np*(di-1);
                bdrylocal = refel.face2local{f};
                % zero those rows in element_matrix
                element_matrix(bdrylocaldof, :) = 0;
                if prob.bc_type{bid}(1) == 'D'
                    for j=1:length(bdrylocaldof)
                        element_matrix(bdrylocaldof(j), bdrylocaldof(j)) = 1;
                    end
                elseif prob.bc_type{bid}(1) == 'N'
                    normal = grid_data.facenormals(:,fid);
                    for j=1:length(bdrylocal)
"*Ddotn_part*"
                        element_matrix(bdrylocaldof(j), ((di-1)*Np+1):(di*Np)) = Ddotn;
                    end
                end
            end
        end
    end
    
    val(st:en) = element_matrix(:);
end
LHS = sparse(I,J,val,dof*dofspernode,dof*dofspernode);";
    
    println(file, content)
    println(file, "");
end

# This is the vector assembly loop.
# It requires code that was generated for this target.
function matlab_linear_file(code, var)
    file = add_generated_file("Linear.m", dir="src");
    # insert the code part into this skeleton
    content = 
"
nnodes= size(grid_data.allnodes,2);
ne  = mesh_data.nel;
Np = refel.Np;
dofspernode = 3;
RHS = zeros(nnodes*dofspernode,1);
glbdof = zeros(Np*dofspernode,1);

%temporary
time = 0;
dt = 1;
var = 0;

% loop over elements
for e=1:ne
    glb = grid_data.loc2glb(:,e)\';
    nodex = grid_data.allnodes(:,glb);
    [detJ, Jac]  = Utils.geometric_factors(refel, nodex);
    wdetj = refel.wg .* detJ;
    
    for di=1:dofspernode
        glbdof(((di-1)*Np + 1):(di*Np)) = (glb - 1).*dofspernode + di;
    end
    
"*code*"
    
    RHS(glbdof) = element_vector;
end
";
    println(file, content);
    
    # boundary condition
    println(file, "
    for i=1:length(grid_data.bdry)
        bdrynode = grid_data.bdry{i};
        for di=1:dofspernode
            bdrydof = (bdrynode-1)*dofspernode + di;
            RHS(bdrydof) = prob.bc_func(di, i);
        end
    end");
    
end

# Nothing here yet. TODO: time stepper
function matlab_stepper_file()
    file = add_generated_file("Stepper.m", dir="src");
end

# Right now this just tries to plot
function matlab_output_file()
    file = add_generated_file("Output.m", dir="src");
    
    if config.dimension == 1
        content = "
gxy = grid_data.allnodes';

X = gxy(:,1);

for figi = 1:total_components
    figure(figi);
    plot(X, result(figi:dofspernode:length(result)));
end
"
    elseif config.dimension == 2
        content = "
gxy = grid_data.allnodes';

X = gxy(:,1);
Y = gxy(:,2);

DT = delaunay(gxy);

for figi = 1:total_components
    figure(figi);
    trisurf(triangulation(DT, X, Y, result(figi:dofspernode:length(result))), 'edgecolor', 'none')
    view(2);
end
"
    else
        content = "% figure output not ready for 3D."
    end
    
    println(file, content)
    println(file, "");
end

###########################################################################################################
# Symbolic to code generation
###########################################################################################################

# If needed, build derivative matrices
function matlabtarget_build_derivative_matrices(lorr, vors)
    code = 
"
% Note on derivative matrices:
% RQn are vandermond matrices for the derivatives of the basis functions
% with Jacobian factors. They are made like this.
% |RQ1|   | rx sx tx || Qx |
% |RQ2| = | ry sy ty || Qy |
% |RQ3|   | rz sz tz || Qz |

"
    if config.dimension == 1
        if vors == "volume"
            code *= "RQ1 = diag(Jac.rx) * refel.Qr;\n";
            code *= "RD1 = diag(Jac.rx) * refel.Ddr;\n";
            code *= "TRQ1 = RQ1';\n"
        else
            # TODO see DG
        end
        
    elseif config.dimension == 2
        if vors == "volume"
            code *= "RQ1 = [diag(Jac.rx) diag(Jac.sx)] * [refel.Qr; refel.Qs];\n";
            code *= "RQ2 = [diag(Jac.ry) diag(Jac.sy)] * [refel.Qr; refel.Qs];\n";
            code *= "RD1 = [diag(Jac.rx) diag(Jac.sx)] * [refel.Ddr; refel.Dds];\n";
            code *= "RD2 = [diag(Jac.ry) diag(Jac.sy)] * [refel.Ddr; refel.Dds];\n";
            code *= "TRQ1 = RQ1';\n"
            code *= "TRQ2 = RQ2';\n"
        else
            # TODO see DG
        end
        
    elseif config.dimension == 3
        if vors == "volume"
            code *= "RQ1 = [diag(Jac.rx) diag(Jac.sx) diag(Jac.tx)] * [refel.Qr; refel.Qs; refel.Qt];\n";
            code *= "RQ2 = [diag(Jac.ry) diag(Jac.sy) diag(Jac.ty)] * [refel.Qr; refel.Qs; refel.Qt];\n";
            code *= "RQ3 = [diag(Jac.rz) diag(Jac.sz) diag(Jac.tz)] * [refel.Qr; refel.Qs; refel.Qt];\n";
            code *= "RD1 = [diag(Jac.rx) diag(Jac.sx) diag(Jac.tx)] * [refel.Ddr; refel.Dds; refel.Ddt];\n";
            code *= "RD2 = [diag(Jac.ry) diag(Jac.sy) diag(Jac.ty)] * [refel.Ddr; refel.Dds; refel.Ddt];\n";
            code *= "RD3 = [diag(Jac.rz) diag(Jac.sz) diag(Jac.tz)] * [refel.Ddr; refel.Dds; refel.Ddt];\n";
            code *= "TRQ1 = RQ1';\n"
            code *= "TRQ2 = RQ2';\n"
            code *= "TRQ3 = RQ3';\n"
        else
            # TODO see DG
        end
    end
    
    return code;
end

# Allocate, compute, or fetch all needed values
function matlabtarget_prepare_needed_values(entities, var, lorr, vors)
    code = "";
    for i=1:length(entities)
        cname = CodeGenerator.make_coef_name(entities[i]);
        if CodeGenerator.is_test_function(entities[i])
            # Assign it a transpose quadrature matrix
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= cname * " = TRQ"*string(entities[i].derivs[di])*"; % d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                end
            else
                code *= cname * " = refel.Q'; % test function.\n";
            end
        elseif CodeGenerator.is_unknown_var(entities[i], var) && lorr == LHS
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= cname * " = RQ"*string(entities[i].derivs[di])*"; % d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                end
            else
                code *= cname * " = refel.Q; % trial function.\n";
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
                    code *= cname * " = 0; % NOTE: derivative applied to constant coefficient = 0\n";
                else
                    code *= cname * " = " * string(cval) * ";\n";
                end
                
            elseif ctype == 2 # a coefficient function
                # This generates something like:
                ######################################
                # coef_n_i = zeros(refel.Np);
                # for coefi = 1:refel.Np
                #     coef_k_i[coefi] = func(x[1,coefi], x[2,coefi],x[3,coefi],time);
                # end
                ######################################
                cargs = "(nodex(coefi), 0, 0, time)";
                if config.dimension == 2
                    cargs = "(nodex(1, coefi), nodex(2, coefi), 0, time)";
                elseif config.dimension == 3
                    cargs = "(nodex(1, coefi), nodex(2, coefi), nodex(3, coefi), time)";
                end
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np,1);\n";cname*" = "*genfunctions[cval].name*"(idx);\n";
                    code *= "for coefi = 1:refel.Np\n\t" * cname * "(coefi) = "*genfunctions[cval].name*"_fun" * cargs * ";\nend\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; % Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        end
                    else
                        code *= cname * " = refel.Q * " * cname * "; % Interpolate at quadrature points.\n";
                    end
                else
                    #TODO surface
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = copy((Finch.variables[1]).values[1, glb])
                if vors == "volume"
                    code *= cname * " = variables("*string(cval)*")("*string(entities[i].index)*", glb);\n";
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; % Apply d/d"*xyzchar[entities[i].derivs[di]]*"\n";
                        end
                    else
                        code *= cname * " = refel.Q * " * cname * "; % Interpolate at quadrature points.\n";
                    end
                else
                    #TODO surface
                end
                
            end
        end # if coefficient
    end # entity loop
    
    return code;
end

function matlabtarget_make_elemental_computation(terms, var, lorr, vors)
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
    
    code = "";
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "element_vector = zeros(refel.Np * "*string(dofsper)*", 1); % Allocate the returned vector.\n"
        else
            code *= "element_matrix = zeros(refel.Np * "*string(dofsper)*", refel.Np * "*string(dofsper)*"); % Allocate the returned matrix.\n"
        end
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
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
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (term_result, test_ind, trial_ind) = matlabtarget_generate_term_calculation(terms[vi][ci][i], var, lorr);
                        
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
                        
                        
                        if length(submatrices[submat_ind]) > 1
                            submatrices[submat_ind] *= " + " * term_result;
                        else
                            submatrices[submat_ind] = term_result;
                        end
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (term_result, test_ind, trial_ind) = matlabtarget_generate_term_calculation(terms[ci][i], var, lorr);
                    
                    # Find the appropriate submatrix for this term
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                    else
                        submat_ind = test_ind;
                    end
                    
                    if length(submatrices[submat_ind]) > 1
                        submatrices[submat_ind] *= " + " * term_result;
                    else
                        submatrices[submat_ind] = term_result;
                    end
                end
            end
            
        end
        
        # Put the submatrices together into element_matrix or element_vector
        if lorr == LHS
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(submatrices[emi, emj]) > 1
                        rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        rangej = "("*string(emj-1)*"*refel.Np + 1):("*string(emj)*"*refel.Np)";
                        code *= "element_matrix("*rangei*", "*rangej*") = " * submatrices[emi,emj] * ";\n";
                    end
                end
            end
            #code *= "return element_matrix;\n"
            
        else # RHS
            for emi=1:dofsper
                if length(submatrices[emi]) > 1
                    rangei = "(("*string(emi-1)*")*refel.Np + 1):("*string(emi)*"*refel.Np)";
                    code *= "element_vector("*rangei*") = " * submatrices[emi] * ";\n";
                end
            end
        end
        #code *= "return element_vector;\n"
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            result = "zeros(refel.Np, refel.Np)";
        else
            result = "zeros(refel.Np,1)";
        end
        
        #process each term
        for i=1:length(terms)
            (term_result, test_ind, trial_ind) = matlabtarget_generate_term_calculation(terms[i], var, lorr);
            
            if i > 1
                result *= " + " * term_result;
            else
                result = term_result;
            end
        end
        if lorr == LHS
            code *= "element_matrix = " * result * ";\n";
        else
            code *= "element_vector = " * result * ";\n";
        end
        #code *= "return " * result * ";\n";
    end
    
    return code;
end

function matlabtarget_generate_term_calculation(term, var, lorr)
    result = "";
    
    if lorr == LHS
        (test_part, trial_part, coef_part, test_ind, trial_ind) = CodeGenerator.separate_factors(term, var);
        # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
        if !(coef_part === nothing)
            result = string(CodeGenerator.replace_entities_with_symbols(test_part)) * " * diag(wdetj .* " * 
                    string(CodeGenerator.replace_entities_with_symbols(coef_part)) * ") * " * 
                    string(CodeGenerator.replace_entities_with_symbols(trial_part));
        else # no coef_part
            result = string(CodeGenerator.replace_entities_with_symbols(test_part)) * " * diag(wdetj) * " * 
                    string(CodeGenerator.replace_entities_with_symbols(trial_part));
        end
    else
        (test_part, trial_part, coef_part, test_ind, trial_ind) = CodeGenerator.separate_factors(term);
        # RHS: test_part * (weight_part .* coef_part)
        if !(coef_part === nothing)
            result = string(CodeGenerator.replace_entities_with_symbols(test_part)) * " * (wdetj .* " * 
                    string(CodeGenerator.replace_entities_with_symbols(coef_part)) * ")";
        else
            result = string(CodeGenerator.replace_entities_with_symbols(test_part)) * " * (wdetj)";
        end
    end
    
    return (result, test_ind, trial_ind);
end
