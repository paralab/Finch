#=
# Format specific functions for writing data output files.
=#

function output_values_raw(vars, file)
    # do one variable/component at a time
    if typeof(vars) <: Array 
        for vi=1:length(vars) # each variable
            output_values_raw(vars[vi], file);
        end
    else
        # Here vars is just one variable
        N = size(vars.values, 2); # number of points
        for ci=1:size(vars.values,1) # each component for that variable
            file.write(vars.values[ci,:]); # write the values
        end
    end
    
end

function output_values_csv(vars, file)
    # columns will be like x, y, z, u1,...
    if typeof(vars) <: Array
        N = size(vars[1].values, 2); # number of points (must be same for all vars)
        if vars[1].location == CELL
            x = fv_info.cellCenters;
        else
            x = grid_data.allnodes;
        end
        
        # header
        line = "x";
        if dim > 1
            line *= ", y";
        end
        if dim > 2
            line *= ", z";
        end
        for vi=1:length(vars) # each variable
            if size(vars[vi].values,1) > 1
                for ci=1:size(vars[vi].values,1)
                    line *= ", "*string(vars[vi].symbol)*"_"*string(ci);
                end
            else
                line *= ", "*string(vars[vi].symbol);
            end
        end
        println(file, line);
        
        # data
        for i=1:N
            # x, y, z
            line = string(x[1,i]);
            if dim > 1
                line *= ", "*string(x[2,i]);
            end
            if dim > 2
                line *= ", "*string(x[3,i]);
            end
            
            # u1, u2, ...
            for vi=1:length(vars) # each variable
                for ci=1:size(vars[vi].values,1)
                    line *= ", "*string(vars[vi].values[ci,i]);
                end
            end
            println(file, line);
        end
        
    else
        N = size(vars.values, 2); # number of points
        if vars.location == CELL
            x = fv_info.cellCenters;
        else
            x = grid_data.allnodes;
        end
        dim = size(x,1);
        
        # header
        line = "x";
        if dim > 1
            line *= ", y";
        end
        if dim > 2
            line *= ", z";
        end
        if size(vars.values,1) > 1
            for ci=1:size(vars.values,1)
                line *= ", "*string(vars.symbol)*"_"*string(ci);
            end
        else
            line *= ", "*string(vars.symbol);
        end
        println(file, line);
        
        # data
        for i=1:N
            # x, y, z
            line = string(x[1,i]);
            if dim > 1
                line *= ", "*string(x[2,i]);
            end
            if dim > 2
                line *= ", "*string(x[3,i]);
            end
            
            # u1, u2, ...
            for ci=1:size(vars.values,1)
                line *= ", "*string(vars.values[ci,i]);
            end
            println(file, line);
        end
    end
end

function output_values_vtk(vars, file, ascii)
    # Use the WriteVTK package.
    # First build MeshCell array
    nel = size(grid_data.loc2glb,2);
    if config.dimension == 1
        np_to_type = Dict([(2,3)]); # line
    elseif config.dimension == 2
        np_to_type = Dict([(3,5),(4,9)]); # triangle, quad
    else
        np_to_type = Dict([(4,10),(8,12)]); # tet, hex
    end
    
    cells = Array{WriteVTK.MeshCell,1}(undef, nel);
    for ei=1:nel
        nvertex = 0;
        for vi=1:length(grid_data.glbvertex[:,ei])
            if grid_data.glbvertex[vi,ei] > 0
                nvertex = nvertex+1;
            else
                break;
            end
        end
        if nvertex < 2
            printerr("Found element with too few vertices when writing VTK output: el="*string(ei)*", vertexmap="*string(grid_data.glbvertex[:,ei]));
            return;
        end
        cell_type = np_to_type[nvertex];
        cells[ei] = WriteVTK.MeshCell(VTKCellType(cell_type), grid_data.glbvertex[:,ei]);
    end
    
    # Initialize the file
    vtkfile = vtk_grid(file, grid_data.allnodes, cells, ascii=ascii);
    
    # Add values
    if typeof(vars) <: Array
        for vi=1:length(vars)
            if length(vars[vi].symvar) > 1
                for ci=1:length(vars[vi].symvar)
                    compname = string(vars[vi].symbol) * "_" * string(ci);
                    vtkfile[compname] = vars[vi].values[ci,:];
                end
            else
                vtkfile[string(vars[vi].symbol)] = vars[vi].values;
            end
        end
    else
        if length(vars.symvar) > 1
            for ci=1:length(vars.symvar)
                compname = string(vars.symbol) * "_" * string(ci);
                vtkfile[compname] = vars.values[ci,:];
            end
        else
            vtkfile[string(vars.symbol)] = vars.values;
        end
    end
    
    return vtk_save(vtkfile);
end
