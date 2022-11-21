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
    if finch_state.config.solver_type == FV
        grid = finch_state.fv_grid;
    else
        grid = finch_state.grid_data;
    end
    # columns will be like x, y, z, u1,...
    if typeof(vars) <: Array
        N = size(vars[1].values, 2); # number of points (must be same for all vars)
        if vars[1].location == CELL
            x = finch_state.fv_info.cellCenters;
        else
            x = grid.allnodes;
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
            x = finch_state.fv_info.cellCenters;
        else
            x = grid.allnodes;
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
    config = finch_state.config;
    if config.solver_type == FV
        grid = finch_state.fv_grid;
        nel = grid.nel_owned; # Only write owned elements
    else
        grid = finch_state.grid_data;
        nel = size(grid.loc2glb,2); # all elements
    end
    # Use the WriteVTK package.
    # First build MeshCell array
    if config.dimension == 1
        np_to_type = Dict([(2,3)]); # line
    elseif config.dimension == 2
        np_to_type = Dict([(3,5),(4,9)]); # triangle, quad
    else
        np_to_type = Dict([(4,10),(8,12)]); # tet, hex
    end
    
    cells = Array{WriteVTK.MeshCell,1}(undef, nel);
    # nodes = fill(-2e20, 3, size(grid.allnodes, 2));
    for ei=1:nel
        nvertex = 0;
        for vi=1:length(grid.glbvertex[:,ei])
            if grid.glbvertex[vi,ei] > 0
                nvertex = nvertex+1;
            else
                break;
            end
        end
        if nvertex < 2
            printerr("Found element with too few vertices when writing VTK output: el="*string(ei)*", vertexmap="*string(grid.glbvertex[:,ei]));
            return;
        end
        cell_type = np_to_type[nvertex];
        cells[ei] = WriteVTK.MeshCell(VTKCellType(cell_type), grid.glbvertex[1:nvertex,ei]);
        # if config.dimension == 1
        #     nodes[1, grid.glbvertex[1:nvertex,ei]] = grid.allnodes[:, grid.glbvertex[1:nvertex,ei]];
        #     nodes[2:3, grid.glbvertex[1:nvertex,ei]] = zeros(2, nvertex);
        # elseif config.dimension == 2
        #     nodes[1:2, grid.glbvertex[1:nvertex,ei]] = grid.allnodes[:, grid.glbvertex[1:nvertex,ei]];
        #     nodes[3, grid.glbvertex[1:nvertex,ei]] = zeros(1, nvertex);
        # else
        #     nodes[:, grid.glbvertex[1:nvertex,ei]] = grid.allnodes[:, grid.glbvertex[1:nvertex,ei]];
        # end
    end
    
    # Let's try all nodes. If it causes problems, uncomment above to remove ghosts
    if config.dimension == 1
        nodes = zeros(3, size(grid.allnodes, 2));
        nodes[1, :] = grid.allnodes[1, :];
    elseif config.dimension == 2
        nodes = zeros(3, size(grid.allnodes, 2));
        nodes[1:2, :] = grid.allnodes[1:2, :];
    else
        nodes = grid.allnodes;
    end
    
    # Initialize the file
    # For multiple partitions write a pvtk file
    if config.num_procs > 1
        # # remove ghost nodes (they should have coordinates -2e20)
        # nextind = 1;
        # for ni=1:size(nodes,2)
        #     if !(nodes[1,ni] == -2e20)
        #         nodes[:,nextind] = nodes[:,ni];
        #         nextind += 1;
        #     end
        # end
        # nodes = nodes[:,1:nextind-1];
        
        vtkfile = pvtk_grid(file, nodes, cells, ascii=ascii, part=config.proc_rank+1, nparts=config.num_procs);
        vtkfile["partition"] = fill(config.partition_index, nel);
    else
        vtkfile = vtk_grid(file, nodes, cells, ascii=ascii);
    end
    
    # Add values
    if typeof(vars) <: Array
        for vi=1:length(vars)
            if vars[vi].location == CELL
                range = 1:grid.nel_owned
            else
                range = 1:size(vars[vi].values,2);
            end
            
            if length(vars[vi].symvar) > 1
                for ci=1:length(vars[vi].symvar)
                    compname = string(vars[vi].symbol) * "_" * string(ci);
                    vtkfile[compname] = vars[vi].values[ci,range];
                end
            else
                vtkfile[string(vars[vi].symbol)] = vars[vi].values[range];
            end
        end
    else
        if vars.location == CELL
            range = 1:grid.nel_owned
        else
            range = 1:size(vars.values,2);
        end
        
        if length(vars.symvar) > 1
            for ci=1:length(vars.symvar)
                compname = string(vars.symbol) * "_" * string(ci);
                vtkfile[compname] = vars.values[ci,range];
            end
        else
            vtkfile[string(vars.symbol)] = vars.values[range];
        end
    end
    
    return vtk_save(vtkfile);
end

########################################################################################################
# My attempt at writing a vtu file.
function output_values_myvtu(var, file, ascii)
    config = finch_state.config;
    if config.solver_type == FV
        grid = finch_state.fv_grid;
    else
        grid = finch_state.grid_data;
    end
    
    if !(typeof(var) <: Array)
        vars = [var];
    else
        vars = var;
    end
    
    cells_only = config.solver_type == FV; # If FV, only use cell data.
    num_points = size(grid.allnodes,2);
    num_cells = grid.nel_owned;
    if ascii
        format = "ascii";
    else
        format = "appended";
    end
    # For appended data
    appended = zeros(UInt64, 0);
    chunk_offset = 1;
    byte_offset = 0;
    
    # cell types: ??
    nodes_per_element = size(grid.glbvertex,1);
    cell_types = zeros(Int8, num_cells);
    cells_per_line = 16;
    if config.dimension == 1
        for ci=1:num_cells
            cell_types[ci] = 3;
        end
    elseif config.dimension == 2
        for ci=1:num_cells
            if nodes_per_element == 3
                cell_types[ci] = 5;
            elseif nodes_per_element == 4
                cell_types[ci] = 9;
            end
        end
    else
        for ci=1:num_cells
            if nodes_per_element == 4
                cell_types[ci] = 10;
            elseif nodes_per_element == 8
                cell_types[ci] = 12;
            end
        end
    end
    
    content = "";
    
    # header
    content *= "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    content *= "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\""*(ascii ? "" : " header_type=\"UInt64\"")*">\n";
    content *= "  <UnstructuredGrid>\n";
    content *= "    <Piece NumberOfPoints=\""*string(num_points)*"\" NumberOfCells=\""*string(num_cells)*"\">\n";
    
    # Points
    content *= "      <Points>\n";
    if ascii
        content *= "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    else
        content *= "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n";
        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_points*3 * 8, UInt64);
    end
    for ni=1:num_points
        scoper = chunk_offset + byte_offset;
        if ascii
            content *= "          ";
            content *= string(grid.allnodes[1, ni]) * " ";
            content *= (config.dimension > 1) ? (string(grid.allnodes[2, ni]) * " ") : "0 ";
            content *= (config.dimension > 2) ? (string(grid.allnodes[3, ni]) * " ") : "0 ";
            content *= "\n";
        else
            tmp1 = grid.allnodes[1, ni];
            tmp2 = (config.dimension > 1) ? grid.allnodes[2, ni] : 0.0;
            tmp3 = (config.dimension > 2) ? grid.allnodes[3, ni] : 0.0;
            (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, tmp1, Float64);
            (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, tmp2, Float64);
            (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, tmp3, Float64);
        end
    end
    content *= "        </DataArray>\n";
    content *= "      </Points>\n";
    
    # Cells
    content *= "      <Cells>\n";
    if ascii
        content *= "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    else
        content *= "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""*string((chunk_offset-1)*8 + byte_offset)*"\">\n";
        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_cells*nodes_per_element * 4, UInt64);
    end
    for ci=1:num_cells
        scoper = chunk_offset + byte_offset;
        if ascii
            content *= "          ";
            for ni=1:nodes_per_element
                content *= string(grid.glbvertex[ni, ci]-1) * " ";
            end
            content *= "\n";
        else
            for ni=1:nodes_per_element
                (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, grid.glbvertex[ni, ci]-1, Int32);
            end
            #println("cell "*string(ci)*" nn "*string(nodes_per_element)*" co "*string(chunk_offset)*" bo "*string(byte_offset))
        end
        
    end
    content *= "        </DataArray>\n";
    
    if ascii
        content *= "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    else
        #println("doing offsets co "*string(chunk_offset)*" bo "*string(byte_offset)*" -> "*string((chunk_offset-1)*8 + byte_offset))
        content *= "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""*string((chunk_offset-1)*8 + byte_offset)*"\">\n";
        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_cells * 4, UInt64);
    end
    tmp = 0;
    offset = 0;
    for ci=1:num_cells
        scoper = chunk_offset + byte_offset;
        if ascii
            if (tmp == 0) content *= "          "; end
            offset += nodes_per_element;
            content *= string(offset) * " ";
            tmp += 1;
            if (tmp == cells_per_line) 
                content *= "\n"; 
                tmp = 0; 
            elseif ci == num_cells
                content *= "\n"; 
            end
        else
            offset += nodes_per_element;
            (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, offset, Int32);
            #println("cell "*string(ci)*" off "*string(offset)*" co "*string(chunk_offset)*" bo "*string(byte_offset))
        end
    end
    content *= "        </DataArray>\n";
    
    if ascii
        content *= "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    else
        #println("doing types co "*string(chunk_offset)*" bo "*string(byte_offset)*" -> "*string((chunk_offset-1)*8 + byte_offset))
        content *= "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""*string((chunk_offset-1)*8 + byte_offset)*"\">\n";
        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_cells, UInt64);
    end
    tmp = 0;
    for ci=1:num_cells
        scoper = chunk_offset + byte_offset;
        if ascii
            if (tmp == 0) content *= "          "; end
            content *= string(cell_types[ci]) * " ";
            tmp += 1;
            if (tmp == cells_per_line) 
                content *= "\n"; 
                tmp = 0; 
            elseif ci == num_cells
                content *= "\n"; 
            end
        else
            (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, cell_types[ci], Int8);
        end
    end
    content *= "        </DataArray>\n";
    content *= "      </Cells>\n";
    
    # Point data
    content *= "      <PointData>\n";
    if !cells_only
        for vi=1:length(vars)
            scoper = chunk_offset + byte_offset;
            comps = size(vars[vi].values,1);
            if ascii
                content *= "        <DataArray type=\"Float64\" Name=\""*string(vars[vi].symbol)*"\""*
                                " NumberOfComponents=\""*string(comps)*"\"" * " format=\"ascii\">\n";
            else
                content *= "        <DataArray type=\"Float64\" Name=\""*string(vars[vi].symbol)*"\""*
                                " NumberOfComponents=\""*string(comps)*"\"" * " format=\"appended\" offset=\""*string((chunk_offset-1)*8 + byte_offset)*"\">\n";
                (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_points*comps * 8, UInt64);
            end
            tmp = 0;
            for ni=1:num_points
                if ascii
                    if (tmp == 0) content *= "          "; end
                    if comps > 1
                        for d=1:comps
                            content *= string(vars[vi].values[d, ni]) * " ";
                        end
                        content *= "\n";
                    else
                        content *= string(vars[vi].values[1, ni]) * " ";
                        tmp += 1;
                        if (tmp >= cells_per_line/2) 
                            content *= "\n"; 
                            tmp = 0;
                        elseif ni==num_points
                            content *= "\n"; 
                        end
                    end
                else
                    for d=1:comps
                        scoper = chunk_offset + byte_offset;
                        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, vars[vi].values[d, ni], Float64);
                    end
                end
            end
            content *= "        </DataArray>\n";
        end
    end
    content *= "      </PointData>\n";
    
    # Cell data
    content *= "      <CellData>\n";
    if cells_only
        for vi=1:length(vars)
            scoper = chunk_offset + byte_offset;
            comps = size(vars[vi].values,1);
            if ascii
                content *= "        <DataArray type=\"Float64\" Name=\""*string(vars[vi].symbol)*
                                " NumberOfComponents=\""*string(comps)*"\"" * " format=\"ascii\">\n";
            else
                content *= "        <DataArray type=\"Float64\" Name=\""*string(vars[vi].symbol)*
                                " NumberOfComponents=\""*string(comps)*"\"" * " format=\"appended\" offset=\""*string((chunk_offset-1)*8 + byte_offset)*"\">\n";
                (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, num_cells*comps * 8, UInt64);
            end
            tmp = 0;
            for ni=1:num_cells
                if ascii
                    if (tmp == 0) content *= "          "; end
                    if comps > 1
                        for d=1:comps
                            content *= string(vars[vi].values[d, ni]) * " ";
                        end
                        content *= "\n";
                    else
                        content *= string(vars[vi].values[1, ni]) * " ";
                        tmp += 1;
                        if (tmp >= cells_per_line/2) 
                            content *= "\n"; 
                            tmp = 0;
                        elseif ni==num_cells
                            content *= "\n"; 
                        end
                    end
                else
                    for d=1:comps
                        scoper = chunk_offset + byte_offset;
                        (chunk_offset, byte_offset) = add_to_appended_data!(appended, chunk_offset, byte_offset, vars[vi].values[d, ni], Float64);
                    end
                end
            end
            content *= "        </DataArray>\n";
        end
    end
    content *= "      </CellData>\n";
    
    content *= "    </Piece>\n";
    content *= "  </UnstructuredGrid>\n";
    
    # appended data if needed
    if !ascii
        content *= "  <AppendedData encoding=\"raw\">\n";
        content *= "_";
        print(file, content);
        for i=1:length(appended)
            write(file, htol(appended[i]));
        end
        
        content = "\n  </AppendedData>\n";
    end
    
    # ending
    content *= "</VTKFile>\n";
    
    println(file, content);
end

# adds one piece of data to an appended data set.
function add_to_appended_data!(data, chunk_offset, byte_offset, newdata, type)
    bytes = UInt64(0);
    tmp = UInt64(0);
    bc = sizeof(type);
    if type == Float64
        intdata = reinterpret(UInt64, newdata)[1];
    elseif type == Float32
        intdata = reinterpret(UInt32, Float32(newdata))[1];
    else
        intdata = newdata;
    end
    
    # insert the bytes of newdata in little endian
    for i=1:bc
        tmp = (intdata >>> (8*(i-1))) & 0xFF;
        tmp = tmp << (8*(i-1));
        bytes = bytes | tmp;
        tmp = UInt64(0);
    end
    
    if byte_offset > 0
        lowbc = min(8-byte_offset, bc);
        highbc = bc-lowbc;
        lowtmp = data[chunk_offset];
        hightmp = UInt64(0);
        
        lowtmp = lowtmp | (bytes << (8*(byte_offset)));
        
        hightmp = bytes >>> (8*(8-byte_offset));
        
        data[chunk_offset] = lowtmp;
        if highbc > 0
            push!(data, hightmp);
            return (chunk_offset + 1, highbc);
        else
            byte_offset = byte_offset + lowbc;
            if byte_offset == 8
                chunk_offset += 1;
                byte_offset = 0;
            end
            return (chunk_offset, byte_offset); 
        end
        
    else
        push!(data, bytes);
        if bc == 8
            chunk_offset += 1;
            byte_offset = 0;
        else
            byte_offset = bc;
        end
            
        return (chunk_offset, byte_offset);
    end
end