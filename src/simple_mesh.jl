# A set of simple mesh makers for testing
export simple_line_mesh, simple_quad_mesh, simple_hex_mesh

#=
# Builds a 1D interval mesh.
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_line_mesh(nx, bn, interval)
    # mesh only
    Nv = nx;                    # number of vertex nodes
    xv = zeros(1,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = nx-1;                 # number of elements
    el = zeros(Int, 2, nel);    # element vertex maps
    etypes = ones(Int, nel);    # element types (gmsh number)
    numvert = 2*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    f2n = zeros(Int, 1,Nv);     # face2node mapping
    f2e = zeros(Int, 2,Nv);     # face2element mapping
    e2f = zeros(Int, 2,nel);    # element2face mapping
    normals = ones(1,Nv);      # normals of faces
    bdryID = zeros(Int, Nv);    # BID of faces (0=interior)
    
    scale = interval[2]-interval[1];
    h = scale/(nx-1);               # uniformly divided
    
    # vertex nodes are initially ordered lexicographically (that's a 17 letter word!)
    for j=1:nx
        xv[1,j] = interval[1] + (j-1)*h;
    end
    
    # Elements are ordered the same way
    for ei=1:nel
        el[1,ei] = ei;
        el[2,ei] = ei+1;
        e2f[1,ei] = ei;
        e2f[2,ei] = ei+1;
    end
    
    # face2node, face2element, normals
    for fi=1:Nv
        f2n[1,fi] = fi;
        f2e[1,fi] = fi-1;
        f2e[2,fi] = fi;
        # normals[fi] = 1; # already set to 1
        # bdryID[fi] = 0; # already set to 0
    end
    f2e[2,Nv] = 0;
    bdryID[1] = 1;
    f2e[:,1] = [1,0];
    normals[1] = -1;
    if bn == 2
        bdryID[Nv] = 2; # for two BIDs
    else
        bdryID[Nv] = 1; # for one BID
    end
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    (refel, grid) = grid_from_mesh(mesh);
    
    # if bn == 2
    #     add_boundary_ID_to_grid(2, x -> x >= interval[2], grid);
    # end
    
    return (mesh, refel, grid);
end

#=
# Builds a 2D quad mesh
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_quad_mesh(nxy, bn, interval)
    if length(nxy) == 2
        nx = nxy[1];
        ny = nxy[2];
    else
        nx = nxy;
        ny = nx;
    end
    if length(interval) == 2
        interval = [interval[1], interval[2], interval[1], interval[2]];
    end
    
    # mesh
    Nv = nx*ny;                 # number of vertex nodes
    xv = zeros(2,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = (nx-1)*(ny-1);        # number of elements
    el = zeros(Int, 4, nel);    # element vertex maps
    etypes = 3*ones(Int, nel);  # element types (gmsh number)
    numvert = 4*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    Nf = nx*(ny-1) + ny*(nx-1); # number of faces
    f2n = zeros(Int, 2,Nf);     # face2node mapping
    f2e = zeros(Int, 2,Nf);     # face2element mapping
    e2f = zeros(Int, 4,nel);    # element2face mapping
    normals = ones(2,Nf);       # normals of faces
    bdryID = zeros(Int, Nf);    # BID of faces (0=interior)
    
    # if bn == 4
    #     bids = [1,2,3,4]; # x=0, x=1, y=0, y=1
    #     allbids = [1,2,3,4];
    # elseif bn == 3
    #     bids = [1,2,3]; # x=0 , x=1, y=0,1
    #     allbids = [1,2,3,3];
    # elseif bn == 2
    #     bids = [1,2]; # x=0,1 , y=0,1
    #     allbids = [1,1,2,2];
    # else
    #     bids = [1]; # everywhere
    #     allbids = [1,1,1,1];
    # end
    bids = [1]; # everywhere
    allbids = [1,1,1,1];
    
    scalex = interval[2]-interval[1];
    hx = scalex/(nx-1);         # uniformly divided
    scaley = interval[4]-interval[3];
    hy = scaley/(ny-1);         # uniformly divided
    
    # vertex nodes are ordered lexicographically
    for j=1:ny
        for i=1:nx
            k = i + (j-1)*nx;
            xv[1, k] = interval[1] + (i-1)*hx;
            xv[2, k] = interval[3] + (j-1)*hy;
        end
    end
    
    # Elements are ordered the same way
    for j=1:(ny-1)
        for i=1:(nx-1)
            ei = i + (j-1)*(nx-1); # element index
            
            el[1,ei] = i + (j-1)*nx;
            el[2,ei] = i + (j-1)*nx + 1;
            el[3,ei] = i + (j)*nx + 1;
            el[4,ei] = i + (j)*nx;
            
            # face2node, face2element, element2face, normals
            f1 = 2*(ei-1)+1 + (j-1); # left face index
            f2 = f1+1; # bottom
            f3 = f1+2; # right
            f4 = j==(ny-1) ? Nf-(nx-1)+i : j*2*(nx-1) + j + i*2; # top
            
            f2n[1,f1] = el[1,ei];
            f2n[2,f1] = el[4,ei];
            f2n[1,f2] = el[1,ei];
            f2n[2,f2] = el[2,ei];
            f2n[1,f3] = el[2,ei];
            f2n[2,f3] = el[3,ei];
            f2n[1,f4] = el[4,ei];
            f2n[2,f4] = el[3,ei];
            
            f2e[2,f1] = ei;
            f2e[2,f2] = ei;
            f2e[1,f3] = ei;
            f2e[1,f4] = ei;
            
            e2f[1,ei] = f1;
            e2f[2,ei] = f2;
            e2f[3,ei] = f3;
            e2f[4,ei] = f4;
            
            normals[:,f1] = [1,0];
            normals[:,f2] = [0,1];
            normals[:,f3] = [1,0];
            normals[:,f4] = [0,1];
        end
    end
    
    # boundaries
    for j=1:(ny-1)
        eleft = (j-1)*(nx-1) + 1;
        eright = j*(nx-1);
        
        fleft = e2f[1,eleft];
        fright = e2f[3,eright];
        
        bdryID[fleft] = 1; # always 1
        bdryID[fright] = allbids[2];
        
        # need to change normals, f2e for left side
        normals[:,fleft] = [-1,0];
        f2e[:,fleft] = [eleft, 0];
    end
    for i=1:(nx-1)
        ebottom = i;
        etop = nel - (nx-1) + i;
        
        fbottom = e2f[2,ebottom];
        ftop = e2f[4,etop];
        
        bdryID[fbottom] = allbids[3];
        bdryID[ftop] = allbids[4];
        
        # need to change normals, f2e for bottom side
        normals[:,fbottom] = [0,-1];
        f2e[:,fbottom] = [ebottom, 0];
    end
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    (refel, grid) = grid_from_mesh(mesh);
    
    if bn == 2
        add_boundary_ID_to_grid(2, (x,y) -> (y <= interval[3]) || (y >= interval[4]), grid);
    elseif bn == 3
        add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid);
        add_boundary_ID_to_grid(3, (x,y) -> ((y <= interval[3]) || (y >= interval[4])) && (x > interval[1] && x < interval[2]), grid);
    elseif bn == 4
        add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid);
        add_boundary_ID_to_grid(3, (x,y) -> (y <= interval[3] && (x > interval[1] && x < interval[2])), grid);
        add_boundary_ID_to_grid(4, (x,y) -> (y >= interval[4] && (x > interval[1] && x < interval[2])), grid);
    end
    
    return (mesh, refel, grid);
end

#=
# Builds a 3D hex mesh
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_hex_mesh(nxyz, bn, interval)
    if length(nxyz) == 3
        nx = nxyz[1];
        ny = nxyz[2];
        nz = nxyz[3]
    else
        nx = nxyz;
        ny = nx;
        nz = nx;
    end
    if length(interval) == 2
        interval = [interval[1], interval[2], interval[1], interval[2], interval[1], interval[2]];
    end
    
    # mesh
    Nv = nx*ny*nz;              # number of vertex nodes
    xv = zeros(3,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = (nx-1)*(ny-1)*(nz-1); # number of elements
    el = zeros(Int, 8, nel);    # element vertex maps
    etypes = 5*ones(Int, nel);  # element types (gmsh number)
    numvert = 8*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    Nf = nz*(nx-1)*(ny-1) + ny*(nx-1)*(nz-1) + nx*(nz-1)*(ny-1); # number of faces
    f2n = zeros(Int, 4,Nf);     # face2node mapping
    f2e = zeros(Int, 2,Nf);     # face2element mapping
    e2f = zeros(Int, 6,nel);    # element2face mapping
    normals = ones(3,Nf);       # normals of faces
    bdryID = zeros(Int, Nf);    # BID of faces (0=interior)
    
    # if bn == 6
    #     bids = [1,2,3,4,5,6]; # all separate
    #     allbids = [1,2,3,4,5,6];
    # elseif bn == 5
    #     bids = [1,2,3,4,5]; # combine z
    #     allbids = [1,2,3,4,5,5];
    # elseif bn == 4
    #     bids = [1,2,3,4]; # combine y and z
    #     allbids = [1,2,3,3,4,4];
    # elseif bn == 3
    #     bids = [1,2,3]; # combine x,y,z
    #     allbids = [1,1,2,2,3,3];
    # elseif bn == 2
    #     bids = [1,2]; # x=0, other
    #     allbids = [1,2,2,2,2,2];
    # else
    #     bids = [1]; # everywhere
    #     allbids = [1,1,1,1,1,1];
    # end
    
    bids = [1]; # everywhere
    allbids = [1,1,1,1,1,1];
    
    scalex = interval[2]-interval[1];
    hx = scalex/(nx-1);         # uniformly divided
    scaley = interval[4]-interval[3];
    hy = scaley/(ny-1);         # uniformly divided
    scalez = interval[6]-interval[5];
    hz = scalez/(nz-1);         # uniformly divided
    
    # vertex nodes are ordered lexicographically
    for k=1:nz
        for j=1:ny
            for i=1:nx
                ni = i + (j-1)*nx + (k-1)*nx*ny;
                xv[1, ni] = interval[1] + (i-1)*hx;
                xv[2, ni] = interval[3] + (j-1)*hy;
                xv[3, ni] = interval[5] + (k-1)*hz;
            end
        end
    end
    
    # Elements are ordered lexicographically
    # Elemental vertices are ordered according to GMSH order
    nextface = 1;
    for k=1:(nz-1)
        for j=1:(ny-1)
            for i=1:(nx-1)
                ei = i + (j-1)*(nx-1) + (k-1)*(ny-1)*(nx-1); # element index
                
                el[1,ei] = i + (j-1)*nx + (k-1)*nx*ny;
                el[2,ei] = i + (j-1)*nx + (k-1)*nx*ny + 1;
                el[3,ei] = i + (j)*nx + (k-1)*nx*ny + 1; #GMSH
                el[4,ei] = i + (j)*nx + (k-1)*nx*ny; #GMSH
                el[5,ei] = i + (j-1)*nx + (k)*nx*ny;
                el[6,ei] = i + (j-1)*nx + (k)*nx*ny + 1;
                el[7,ei] = i + (j)*nx + (k)*nx*ny + 1; #GMSH
                el[8,ei] = i + (j)*nx + (k)*nx*ny; #GMSH
                
                # face2node, face2element, element2face, normals
                if i==1
                    f1 = nextface;
                    nextface = nextface+1;
                else
                    f1 = e2f[2,ei-1];# left face index
                end
                f2 = nextface; # right
                nextface = nextface+1;
                
                if j==1
                    f3 = nextface;
                    nextface = nextface+1;
                else
                    f3 = e2f[4,ei-(nx-1)]; # bottom
                end
                f4 = nextface; # top
                nextface = nextface+1;
                
                if k==1
                    f5 = nextface;
                    nextface = nextface+1;
                else
                    f5 = e2f[6,ei-((nx-1)*(ny-1))]; # front
                end
                f6 = nextface; # back
                nextface = nextface+1;
                
                f2n[:,f1] = el[[1,4,5,8],ei]; # modified to match GMSH
                f2n[:,f2] = el[[2,3,6,7],ei];
                f2n[:,f3] = el[[1,2,5,6],ei];
                f2n[:,f4] = el[[4,3,7,8],ei];
                f2n[:,f5] = el[[1,2,3,4],ei];
                f2n[:,f6] = el[[5,6,7,8],ei];
                
                f2e[2,f1] = ei;
                f2e[1,f2] = ei;
                f2e[2,f3] = ei;
                f2e[1,f4] = ei;
                f2e[2,f5] = ei;
                f2e[1,f6] = ei;
                
                e2f[1,ei] = f1;
                e2f[2,ei] = f2;
                e2f[3,ei] = f3;
                e2f[4,ei] = f4;
                e2f[5,ei] = f5;
                e2f[6,ei] = f6;
                
                normals[:,f1] = [1,0,0];
                normals[:,f2] = [1,0,0];
                normals[:,f3] = [0,1,0];
                normals[:,f4] = [0,1,0];
                normals[:,f5] = [0,0,1];
                normals[:,f6] = [0,0,1];
            end
        end
    end
    
    # boundaries
    for k=1:(nz-1)
        for j=1:(ny-1)
            eleft = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
            eright = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
            
            fleft = e2f[1,eleft];
            fright = e2f[2,eright];
            
            bdryID[fleft] = 1; # always 1
            bdryID[fright] = allbids[2];
            
            # need to change normals, f2e for left side
            normals[:,fleft] = [-1,0,0];
            f2e[:,fleft] = [eleft, 0];
        end
    end
    
    for k=1:(nz-1)
        for i=1:(nx-1)
            ebottom = (k-1)*(nx-1)*(ny-1) + i;
            etop = (k-1)*(nx-1)*(ny-1) + (nx-1)*(ny-2) + i;
            
            fbottom = e2f[3,ebottom];
            ftop = e2f[4,etop];
            
            bdryID[fbottom] = allbids[3];
            bdryID[ftop] = allbids[4];
            
            # need to change normals, f2e for bottom side
            normals[:,fbottom] = [0,-1, 0];
            f2e[:,fbottom] = [ebottom, 0];
        end
    end
    
    for j=1:(ny-1)
        for i=1:(nx-1)
            efront = (j-1)*(nx-1) + i;
            eback = (nz-2)*(nx-1)*(ny-1) + efront;
            
            ffront = e2f[5,efront];
            fback = e2f[6,eback];
            
            bdryID[ffront] = allbids[5];
            bdryID[fback] = allbids[6];
            
            # need to change normals, f2e for bottom side
            normals[:,ffront] = [0,0,-1];
            f2e[:,ffront] = [efront, 0];
        end
    end
    
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    (refel, grid) = grid_from_mesh(mesh);
    
    tiny = 1e-13;
    if bn == 6
        # bids = [1,2,3,4,5,6]; # all separate
        add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid);
        add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid);
        add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid);
        add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny), grid);
        add_boundary_ID_to_grid(6, (x,y,z) -> (z >= interval[6]-tiny), grid);
    elseif bn == 5
        # bids = [1,2,3,4,5]; # combine z
        add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid);
        add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid);
        add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid);
        add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid);
    elseif bn == 4
        # bids = [1,2,3,4]; # combine y and z
        add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]), grid);
        add_boundary_ID_to_grid(3, (x,y,z) -> ((y <= interval[3]+tiny) || (y >= interval[4]-tiny)), grid);
        add_boundary_ID_to_grid(4, (x,y,z) -> ((z <= interval[5]+tiny) || (z >= interval[6]-tiny)), grid);
    elseif bn == 3
        # bids = [1,2,3]; # combine x,y,z
        add_boundary_ID_to_grid(2, (x,y,z) -> (y <= interval[3]+tiny) || (y >= interval[4]-tiny), grid);
        add_boundary_ID_to_grid(3, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid);
    elseif bn == 2
        # bids = [1,2]; # x=0, other
        add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny) || (y <= interval[3]+tiny) || (y >= interval[4]-tiny) || (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid);
    end
    
    return (mesh, refel, grid)
end
