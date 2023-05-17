#=
Makes a mixed mesh with quads and triangles.
Writes it to a Medit file
=#

## Parameters ##############################################
# Depth of coarsest level
min_depth = 4;

# Refinement criterion
# An element centered at (x,y) is refined if true
function refine_region(x, y)
    r = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
    return (x <= 0.05) || (x >= 0.95) || (y <= 0.05) || (y >= 0.95) || (r < 0.2)
    # return x < 0.05 || y < 0.05 || y > 0.95 || x > 0.6
end
############################################################

#= The configuration of neighbors c = coarse, f = fine
  c           c           c           c           f
c   c = 0   f   c = 1   c   c = 2   c   f = 3   c   c = 4
  c           c           f           c           c
  
  c           c           f           f
f   c = 5   c   f = 6   c   f = 7   f   c = 8
  f           f           c           c
  
  c           f
f   f = 9   c   c = 10
  c           f
  
  c           f           f           f
f   f = 11  c   f = 12  f   f = 13  f   c = 14
  f           f           c           f

  f
f   f = 15
  f

Really, 9 to 15 should not happen.
=#
function get_config(to_refine, neighbors, eid)
    if to_refine[eid]
        return 16;
    end
    map = [1,2,4,8, 3,6,12,9, 5,10, 11,7,14,13, 15];
    imap = zeros(Int,15);
    for i=1:15
        imap[map[i]] = i;
    end
    
    config = 0;
    if neighbors[1,eid]>0 && to_refine[neighbors[1,eid]]
        config += 1;
    end
    if neighbors[2,eid]>0 && to_refine[neighbors[2,eid]]
        config += 2;
    end
    if neighbors[3,eid]>0 && to_refine[neighbors[3,eid]]
        config += 4;
    end
    if neighbors[4,eid]>0 && to_refine[neighbors[4,eid]]
        config += 8;
    end
    
    if config > 0
        config = imap[config];
    end
    
    return config;
end

# Builds the mesh.
function build_mesh(min_depth)
    # make base mesh
    N1d = 2^min_depth;
    nel_base = N1d * N1d;
    nnodes_base = (N1d+1)*(N1d+1);
    xy_base = zeros(2, nnodes_base);
    elements_base = zeros(Int, 4, nel_base);
    to_refine = zeros(Bool, nel_base);
    neighbors = zeros(Int, 4, nel_base);
    
    # node coordinates
    next = 1;
    for yi=1:(N1d+1)
        for xi=1:(N1d+1)
            xy_base[1, next] = (xi-1) / N1d;
            xy_base[2, next] = (yi-1) / N1d;
            next += 1;
        end
    end
    
    # elements
    next = 1;
    for yi=1:N1d
        for xi=1:N1d
            xoffset = (yi-1)*(N1d+1);
            elements_base[1, next] = xoffset + xi;
            elements_base[2, next] = xoffset + xi + 1;
            elements_base[3, next] = xoffset + xi + N1d + 2;
            elements_base[4, next] = xoffset + xi + N1d + 1;
            
            neighbors[1, next] = (xi>1) ? next-1 : 0;
            neighbors[2, next] = (yi>1) ? next-N1d : 0;
            neighbors[3, next] = (xi<N1d) ? next+1 : 0;
            neighbors[4, next] = (yi<N1d) ? next+N1d : 0;
            next += 1;
        end
    end
    
    # return (nel_base, xy_base, elements_base);
    
    # flag for refinement in fine region
    for ei=1:nel_base
        xc = 0.0;
        yc = 0.0;
        for ni=1:4
            nid = elements_base[ni, ei];
            xc = xc + xy_base[1, nid];
            yc = yc + xy_base[2, nid];
        end
        xc /= 4;
        yc /= 4;
        
        to_refine[ei] = refine_region(xc, yc)
    end
    
    # count new nel and nnodes
    nel = nel_base;
    nnodes = nnodes_base;
    for ei=1:nel_base
        # if this is refined, add 3
        if to_refine[ei]
            nel += 3;
            nnodes += 5; # there are duplicates
        else
            # if this has refined neighbors, add as needed
            neigh_config = get_config(to_refine, neighbors, ei);
            
            if neigh_config == 0
                # nothing
            elseif neigh_config < 5
                nel += 2;
            elseif neigh_config < 9
                nel += 3;
            elseif neigh_config < 11
                nel += 1;
            elseif neigh_config < 15
                nel += 3; 
            end
        end
    end
    
    # do the refinement
    xy = zeros(2,nnodes);
    elements = zeros(Int, 4, nnodes);
    new_nodes = zeros(Int, 5, nel_base);
    elements[:, 1:nel_base] .= elements_base;
    xy[:, 1:nnodes_base] .= xy_base;
    
    eltype = fill(3, nel);
    nextel = nel_base+1;
    nextnode = nnodes_base+1;
    # Do one pass refining quads and creating all new nodes
    for yi=1:N1d
        for xi=1:N1d
            eid = (yi-1)*N1d + xi;
            if to_refine[eid]
                if neighbors[1,eid] > 0 && to_refine[neighbors[1,eid]]
                    new_nodes[1,eid] = new_nodes[3,neighbors[1,eid]];
                else
                    new_nodes[1,eid] = nextnode;
                    xy[1,nextnode] = xy[1,elements[1,eid]];
                    xy[2,nextnode] = 0.5*(xy[2,elements[1,eid]] + xy[2,elements[4,eid]]);
                    nextnode += 1;
                end
                if neighbors[2,eid] > 0 && to_refine[neighbors[2,eid]]
                    new_nodes[2,eid] = new_nodes[4,neighbors[2,eid]];
                else
                    new_nodes[2,eid] = nextnode;
                    xy[2,nextnode] = xy[2,elements[1,eid]];
                    xy[1,nextnode] = 0.5*(xy[1,elements[1,eid]] + xy[1,elements[2,eid]]);
                    nextnode += 1;
                end
                new_nodes[3,eid] = nextnode;
                xy[1,nextnode] = xy[1,elements[2,eid]];
                xy[2,nextnode] = 0.5*(xy[2,elements[2,eid]] + xy[2,elements[3,eid]]);
                new_nodes[4,eid] = nextnode+1;
                xy[1,nextnode+1] = 0.5*(xy[1,elements[1,eid]] + xy[1,elements[2,eid]]);
                xy[2,nextnode+1] = xy[2,elements[3,eid]];
                new_nodes[5,eid] = nextnode+2;
                xy[1,nextnode+2] = 0.5*(xy[1,elements[1,eid]] + xy[1,elements[2,eid]]);
                xy[2,nextnode+2] = 0.5*(xy[2,elements[2,eid]] + xy[2,elements[3,eid]]);
                nextnode += 3;
                
                elements[1, nextel] = new_nodes[2,eid];
                elements[2, nextel] = elements[2,eid];
                elements[3, nextel] = new_nodes[3,eid];
                elements[4, nextel] = new_nodes[5,eid];
                eltype[nextel] = 3;
                nextel += 1;
                
                elements[1, nextel] = new_nodes[5,eid];
                elements[2, nextel] = new_nodes[3,eid];
                elements[3, nextel] = elements[3,eid];
                elements[4, nextel] = new_nodes[4,eid];
                eltype[nextel] = 3;
                nextel += 1;
                
                elements[1, nextel] = new_nodes[1,eid];
                elements[2, nextel] = new_nodes[5,eid];
                elements[3, nextel] = new_nodes[4,eid];
                elements[4, nextel] = elements[4,eid];
                eltype[nextel] = 3;
                nextel += 1;
                
                # elements[1, eid] = elements[1,eid];
                elements[2, eid] = new_nodes[2,eid];
                elements[3, eid] = new_nodes[5,eid];
                elements[4, eid] = new_nodes[1,eid];
                eltype[eid] = 3;
            end
        end
    end
    
    # Do a second pass refining triangles (no new nodes)
    for yi=1:N1d
        for xi=1:N1d
            eid = (yi-1)*N1d + xi;
            if !to_refine[eid]
                n_config = get_config(to_refine, neighbors, eid);
                
                if n_config == 0
                    continue;
                elseif n_config == 1
                    new_node = new_nodes[3,neighbors[1,eid]];
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = elements[3,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[3, eid] = new_node;
                    eltype[eid] = 2;
                    
                elseif n_config == 2
                    new_node = new_nodes[4,neighbors[2,eid]];
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = elements[3,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[2, eid] = new_node;
                    elements[3, eid] = elements[4,eid];
                    eltype[eid] = 2;
                    
                elseif n_config == 3
                    new_node = new_nodes[1,neighbors[3,eid]];
                    elements[1, nextel] = elements[1,eid];
                    elements[2, nextel] = elements[2,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[2, eid] = elements[4,eid];
                    elements[3, eid] = new_node;
                    eltype[eid] = 2;
                    
                elseif n_config == 4
                    new_node = new_nodes[2,neighbors[4,eid]];
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = elements[3,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[4,eid];
                    elements[2, nextel] = elements[1,eid];
                    elements[3, nextel] = new_node;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[3, eid] = new_node;
                    eltype[eid] = 2;
                
                elseif n_config == 5
                    new_node1 = new_nodes[3,neighbors[1,eid]];
                    new_node2 = new_nodes[4,neighbors[2,eid]];
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = elements[3,eid];
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node1;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = new_node1;
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[2, eid] = new_node2;
                    elements[3, eid] = new_node1;
                    eltype[eid] = 2;
                    
                elseif n_config == 6
                    new_node1 = new_nodes[4,neighbors[2,eid]];
                    new_node2 = new_nodes[1,neighbors[3,eid]];
                    elements[1, nextel] = elements[3,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[1,eid];
                    elements[2, nextel] = elements[4,eid];
                    elements[3, nextel] = new_node1;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[4,eid];
                    elements[2, nextel] = new_node1;
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, eid] = new_node2;
                    elements[3, eid] = new_node1;
                    eltype[eid] = 2;
                    
                elseif n_config == 7
                    new_node1 = new_nodes[1,neighbors[3,eid]];
                    new_node2 = new_nodes[2,neighbors[4,eid]];
                    elements[1, nextel] = elements[4,eid];
                    elements[2, nextel] = elements[1,eid];
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[1,eid];
                    elements[2, nextel] = elements[2,eid];
                    elements[3, nextel] = new_node1;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[1,eid];
                    elements[2, nextel] = new_node1;
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, eid] = new_node2;
                    elements[2, eid] = new_node1;
                    eltype[eid] = 2;
                    
                elseif n_config == 8
                    new_node1 = new_nodes[2,neighbors[4,eid]];
                    new_node2 = new_nodes[3,neighbors[1,eid]];
                    elements[1, nextel] = elements[1,eid];
                    elements[2, nextel] = elements[2,eid];
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = elements[3,eid];
                    elements[3, nextel] = new_node1;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, nextel] = elements[2,eid];
                    elements[2, nextel] = new_node1;
                    elements[3, nextel] = new_node2;
                    eltype[nextel] = 2;
                    nextel += 1;
                    
                    elements[1, eid] = elements[4, eid]
                    elements[2, eid] = new_node2;
                    elements[3, eid] = new_node1;
                    eltype[eid] = 2;
                
                end
            end
        end
    end
    
    # Cut off the excess
    elements = elements[:, 1:(nextel-1)];
    xy = xy[:, 1:(nextnode-1)];
    
    return (nel, xy, elements, eltype);
end

function mesh_write_medit(file, nel, nodes, elements, etype)
    dim = 2; # size(nodes,1);
    nnodes = size(nodes,2);
    # etype = 3; # quads
    # nv = 4;
    
    println(file, "MeshVersionFormatted 1");
    println(file, "# Created by mesh_read_bte");
    println(file, "");
    println(file, "Dimension");
    println(file, string(dim));
    
    println(file, "Vertices");
    println(file, string(nnodes));
    for i=1:nnodes
        for j=1:dim
            print(file, string(nodes[j,i]) * "    ");
        end
        println(file, string(i));
    end
    
    # Count triangles and quads
    tri_count = 0;
    quad_count = 0;
    for i=1:nel
        if etype[i] == 2
            tri_count += 1;
        else
            quad_count += 1;
        end
    end
    
    # First do quads
    if quad_count > 0
        println(file, "Quadrilaterals");
        println(file, string(quad_count));
        for i=1:nel
            if etype[i] == 3
                for j=1:4
                    print(file, string(elements[j,i]) * "    ");
                end
                println(file, "0");
            end
        end
    end
    
    # Then do triangles
    if tri_count > 0
        println(file, "Triangles");
        println(file, string(tri_count));
        for i=1:nel
            if etype[i] == 2
                for j=1:3
                    print(file, string(elements[j,i]) * "    ");
                end
                println(file, "0");
            end
        end
    end
    
    println(file, "End");
end

# build mesh
(nel, nodes, elements, etype) = build_mesh(min_depth);

println("Built mixed mesh with $nel elements.")

# Write to file
filename = "mixedmesh.mesh"
outfile = open(filename, "w");
mesh_write_medit(outfile, nel, nodes, elements, etype);
close(outfile);

println("Wrote mesh to $filename in MEDIT format");

# Plot to check
println("Plotting... This could take a minute the first time.");
try
    using Plots
    pyplot();
catch e
    println("You don't have the Plots package yet. I'll install it for you now.")
    using Pkg
    Pkg.add("Plots")
    Pkg.add("PyPlot")
    using Plots
    pyplot();
end
elx = zeros(4, nel);
ely = zeros(4, nel);
centers = zeros(2, nel);
for i=1:nel
    if etype[i] == 3
        elx[:,i] = nodes[1, elements[1:4,i]];
        ely[:,i] = nodes[2, elements[1:4,i]];
        centers[1,i] = sum(elx[:,i])/4;
        centers[2,i] = sum(ely[:,i])/4;
    else
        elx[1:3,i] = nodes[1, elements[1:3,i]];
        ely[1:3,i] = nodes[2, elements[1:3,i]];
        elx[4,i] = nodes[1, elements[1,i]];
        ely[4,i] = nodes[2, elements[1,i]];
        centers[1,i] = sum(elx[1:3,i])/3;
        centers[2,i] = sum(ely[1:3,i])/3;
    end
end
p1 = plot(elx, ely, legend=false)
# for i=1:nel
#     annotate!(centers[1,i], centers[2,i], string(i));
# end
display(plot(p1))