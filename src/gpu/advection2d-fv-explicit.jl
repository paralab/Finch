# This is the first GPU supported Finch code; this helps us to create the GPU target for Finch.
# This version is 10x faster compared to the CPU version of Finch.

#=

using Finch
initFinch("FVadvection2d");
useLog("FVadvection2dlog", level=3)

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_IMPLICIT)
setSteps(0.001, 50)

#mesh
use_unstructured=false;
mesh(QUADMESH, elsperdim=[100, 100], bids=4, interval=[0, 0.1, 0, 0.3], partitions=0)

# Variables and BCs
u = variable("u", location=CELL)
boundary(u, 1, FLUX, "(abs(y-0.06) < 0.033 && sin(3*pi*t)>0) ? 1 : 0") # x=0
boundary(u, 2, NO_BC) # x=0.1
boundary(u, 3, NO_BC) # y=0
boundary(u, 4, NO_BC) # y=0.3

initial(u, "0")

# Coefficients
coefficient("a", ["0.1*cos(pi*x/2/0.1)","0.3*sin(pi*x/2/0.1)"], type=VECTOR) # advection velocity
coefficient("s", ["0.1 * sin(pi*x)^4 * sin(pi*y)^4"]) # source

conservationForm(u, "s + surface(upwind(a,u))");
solve(u)
=#

using CUDA

function generated_solve_function_for_u_cuda(var::Vector{Variable{FT}}, mesh::Grid, refel::Refel, geometric_factors::GeometricFactors, fv_info::FVInfo, config::FinchConfig, coefficients::Vector{Coefficient}, variables::Vector{Variable{FT}}, test_functions::Vector{Coefficient}, ordered_indexers::Vector{Indexer}, prob::FinchProblem, time_stepper::Stepper, buffers::ParallelBuffers, timer_output::TimerOutput, nl_var=nothing) where FT<:AbstractFloat

    dofs_per_node = 1;
    num_elements = mesh.nel_owned;
    num_elements_ghost = mesh.nel_ghost;
    num_faces = mesh.nface_owned + mesh.nface_ghost;
    fv_dofs_partition = dofs_per_node * (num_elements + num_elements_ghost);

    c_cellCenter = CuArray(fv_info.cellCenters)
    c_face2element = CuArray(mesh.face2element)
    c_element2face = CuArray(mesh.element2face)
    c_faceCenters = CuArray(fv_info.faceCenters)
    c_faceNormals = CuArray(mesh.facenormals)
    c_geometric_factors_vol = CuArray(geometric_factors.volume)
    c_geometric_factors_area = CuArray(geometric_factors.area)

    bc_type = zeros(Int64, num_faces)
    _fbid = zeros(Int64, num_faces)
    for fid=1:num_faces
        fbid = mesh.facebid[fid]
        if fbid > 0
            _fbid[fid] = 1
            bc_type[fid] = prob.bc_type[1, fbid] == FLUX ? 1 : 0
        end
    end

    c_values = CUDA.zeros(Float64, num_elements)
    c_index = 1
    c_bc_type = CuArray(bc_type)


    c_global_vector = CUDA.zeros(Float64, num_elements)
    t = 0.0
    dt = time_stepper.dt
    c_dt = CUDA.fill(Float64(dt), num_elements)

    for ti = 1:time_stepper.Nsteps #time_stepper.Nsteps
        @cuda threads=256 blocks=ceil(Int, num_elements/256) kernel_elm_loop!(c_global_vector, c_element2face, 
                                                                            c_face2element, c_geometric_factors_area, 
                                                                            c_geometric_factors_vol, c_values, c_index, c_faceCenters, 
                                                                            c_faceNormals, c_cellCenter, c_bc_type, t, num_elements) 
        synchronize()
        @cuda threads=256 blocks=ceil(Int, num_elements/256) kernel_value_update!(c_values, c_global_vector, dt, num_elements)
        synchronize()
        t = (t + dt)
    end
    solution = Array(c_values);
    var[1].values[:] .= solution[:]
    
end # function

function kernel_value_update!(c_values, c_global_vector, c_dt, num_elements)
    eid = threadIdx().x + (blockIdx().x-1) * blockDim().x;
    if eid <= num_elements
        c_values[eid] = c_values[eid]  + c_dt * c_global_vector[eid]
    end
    return nothing
end

function kernel_elm_loop!(c_global_vector, c_element2face, c_face2element, c_geometric_factors_area, c_geometric_factors_vol, 
                            c_values, c_index, c_faceCenters, c_faceNormals, c_cellCenter, c_bc_type, t, num_elements)  ### 1:num_elm
    eid = threadIdx().x + (blockIdx().x-1) * blockDim().x;
    if eid <= num_elements
        # source value
        x = c_cellCenter[1,eid]
        y = c_cellCenter[2,eid]
        value__s_1 = (0.1 * sin(pi*x)^4 * sin(pi*y)^4)
        
        fid1 = c_element2face[1, eid]  # element2face mesh
        fid2 = c_element2face[2, eid]
        fid3 = c_element2face[3, eid]
        fid4 = c_element2face[4, eid]
        
        # a dot n
        xs = 0.1*cos(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[1, fid1]
        ys = 0.3*sin(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[2, fid1]
        adotn1 = xs + ys
        
        xs = 0.1*cos(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[1, fid2]
        ys = 0.3*sin(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[2, fid2]
        adotn2 = xs + ys
        
        xs = 0.1*cos(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[1, fid3]
        ys = 0.3*sin(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[2, fid3]
        adotn3 = xs + ys
        
        xs = 0.1*cos(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[1, fid4]
        ys = 0.3*sin(pi*c_faceCenters[1, fid1]/2/0.1)*c_faceNormals[2, fid4]
        adotn4 = xs + ys
        
        volume = c_geometric_factors_vol[eid] 
        area_over_volume_1 = c_geometric_factors_area[fid1]/volume
        area_over_volume_2 = c_geometric_factors_area[fid2]/volume
        area_over_volume_3 = c_geometric_factors_area[fid3]/volume
        area_over_volume_4 = c_geometric_factors_area[fid4]/volume
        
        left_el1 =  c_face2element[1, fid1]  # face2element  mesh
        right_el1 = c_face2element[2, fid1]
        left_el2 =  c_face2element[1, fid2]
        right_el2 = c_face2element[2, fid2]
        left_el3 =  c_face2element[1, fid3]
        right_el3 = c_face2element[2, fid3]
        left_el4 =  c_face2element[1, fid4]
        right_el4 = c_face2element[2, fid4]
        
        # check for boundary faces
        is_bdry1 = right_el1==0 ? 1 : 0
        is_bdry2 = right_el2==0 ? 1 : 0
        is_bdry3 = right_el3==0 ? 1 : 0
        is_bdry4 = right_el4==0 ? 1 : 0
        
        # switch normal direction if needed
        adotn1 = eid==left_el1 ? adotn1 : -adotn1
        adotn2 = eid==left_el2 ? adotn2 : -adotn2
        adotn3 = eid==left_el3 ? adotn3 : -adotn3
        adotn4 = eid==left_el4 ? adotn4 : -adotn4
        
        # if right is 0, make it equal to left
        right_el1 = right_el1==0 ? left_el1 : right_el1
        right_el2 = right_el2==0 ? left_el2 : right_el2
        right_el3 = right_el3==0 ? left_el3 : right_el3
        right_el4 = right_el4==0 ? left_el4 : right_el4
        
        neighbor1 = eid == right_el1 ? left_el1 : right_el1
        neighbor2 = eid == right_el2 ? left_el2 : right_el2
        neighbor3 = eid == right_el3 ? left_el3 : right_el3
        neighbor4 = eid == right_el4 ? left_el4 : right_el4
        
        value_CELL1_u_1 =   c_values[eid] # variables
        value_CELL2_u_1_1 = c_values[neighbor1]
        value_CELL2_u_1_2 = c_values[neighbor2]
        value_CELL2_u_1_3 = c_values[neighbor3]
        value_CELL2_u_1_4 = c_values[neighbor4]
        
        upwind1 = adotn1>=0 ? value_CELL1_u_1 : value_CELL2_u_1_1
        upwind2 = adotn2>=0 ? value_CELL1_u_1 : value_CELL2_u_1_2
        upwind3 = adotn3>=0 ? value_CELL1_u_1 : value_CELL2_u_1_3
        upwind4 = adotn4>=0 ? value_CELL1_u_1 : value_CELL2_u_1_4
        
        upwind1 = upwind1 * adotn1;
        upwind2 = upwind2 * adotn2;
        upwind3 = upwind3 * adotn3;
        upwind4 = upwind4 * adotn4;

        y1 = c_faceCenters[2, fid1]  # faceCenters fv_info
        y2 = c_faceCenters[2, fid2]
        y3 = c_faceCenters[2, fid3]
        y4 = c_faceCenters[2, fid4]

        flux_case1 = c_bc_type[fid1]
        flux_case2 = c_bc_type[fid2]
        flux_case3 = c_bc_type[fid3]
        flux_case4 = c_bc_type[fid4]

        bdry_flux1 = (abs(y1-0.06) < 0.033 && sin(3*pi*t)>0) ? 1.0 : 0.0
        bdry_flux2 = (abs(y2-0.06) < 0.033 && sin(3*pi*t)>0) ? 1.0 : 0.0
        bdry_flux3 = (abs(y3-0.06) < 0.033 && sin(3*pi*t)>0) ? 1.0 : 0.0
        bdry_flux4 = (abs(y4-0.06) < 0.033 && sin(3*pi*t)>0) ? 1.0 : 0.0
        
        bdry_flux1 = bdry_flux1 * is_bdry1 * flux_case1
        bdry_flux2 = bdry_flux2 * is_bdry2 * flux_case2
        bdry_flux3 = bdry_flux3 * is_bdry3 * flux_case3
        bdry_flux4 = bdry_flux4 * is_bdry4 * flux_case4
        
        _flux1 = (bdry_flux1 - upwind1 * (1 - is_bdry1 * flux_case1)) * area_over_volume_1  # area_over_volume
        _flux2 = (bdry_flux2 - upwind2 * (1 - is_bdry2 * flux_case2)) * area_over_volume_2
        _flux3 = (bdry_flux3 - upwind3 * (1 - is_bdry3 * flux_case3)) * area_over_volume_3
        _flux4 = (bdry_flux4 - upwind4 * (1 - is_bdry4 * flux_case4)) * area_over_volume_4
        
        c_global_vector[eid] = _flux1 + _flux2 + _flux3 + _flux4 + value__s_1
        
    end
    return nothing
end
