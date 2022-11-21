#=
# Create a recursive ordering
# Each state maps to a configuration of states with a rule and ordering
#
# For example, 2D Hilbert has 4 states: n, ], [, u labeled here as s1,s2,s3,s4
#
#       s1 s1
# s1 -> s2 s3  This state, s1, corresponds to rule [s2,s3,s1,s1] and ordering [1,3,4,2]
# 
# Morton is the simplest with one rule [s1,s2,s3,s4] and ordering [1,2,3,4]
=#
export get_recursive_order, reorder_grid_recursive!

struct RecursiveOrdering
    states::Int                 # Number of possible states
    rules::Array{Array{Int,1}}  # Rules for each state mapping sectors to states
    orders::Array{Array{Int,1}} # Ordering for each state
end

# Use this function to build an ordering for a specified grid size
function get_recursive_order(type, dim, griddim)
    nnodes = 1;
    for i=1:length(griddim)
        nnodes = nnodes*griddim[i];
    end
    
    if dim == 1
        # This doesn't make sense, to just keep regular ordering
        return 1:nnodes;
        
    elseif dim == 2
        # Find the smallest N such that 2^N >= max(griddim)
        maxgriddim = max(griddim[1],griddim[2]);
        N=0;
        while maxgriddim > (2^N)
            N = N+1;
        end
        twoN = 2^N;
        
        if type == "hilbert"
            rorder_rules = hilbert_ordering_rules(2); # The RecursiveOrdering object
            rorder = build_ordering_2d(rorder_rules, N, false); # The ordering (false -> inverse)
        elseif type == "morton"
            rorder_rules = morton_ordering_rules(2); # The RecursiveOrdering object
            rorder = build_ordering_2d(rorder_rules, N, false); # The ordering (false -> inverse)
        end
        
        # Truncate ordering to match griddim size
        subrorder = zeros(Int,nnodes);
        subind = 1;
        for i=1:length(rorder)
            # Get coordinates in 2^N grid
            hind = rorder[i];
            hx = mod(hind-1,twoN) + 1;
            hy = Int(floor((hind-1) / (twoN))) + 1;
            
            if hx <= griddim[1] && hy <= griddim[2]
                gind = hx + griddim[1]*(hy-1);
                subrorder[subind] = gind;
                subind += 1;
            end
        end
        return subrorder;
        
    elseif dim == 3
        # Find the smallest N such that 2^N >= max(griddim)
        maxgriddim = max(griddim[1],max(griddim[2],griddim[3]));
        N=0;
        while maxgriddim > (2^N)
            N = N+1;
        end
        twoN = 2^N;
        
        if type == HILBERT_ORDERING
            rorder_rules = hilbert_ordering_rules(3); # The RecursiveOrdering object
            rorder = build_ordering_3d(rorder_rules, N, false); # The ordering (false -> inverse)
        elseif type == MORTON_ORDERING
            rorder_rules = morton_ordering_rules(3); # The RecursiveOrdering object
            rorder = build_ordering_3d(rorder_rules, N, false); # The ordering (false -> inverse)
        end
        
        # Truncate ordering to match griddim size
        subrorder = zeros(Int,nnodes);
        subind = 1;
        for i=1:length(rorder)
            # Get coordinates in 2^N grid
            hind = rorder[i];
            hx = mod(hind-1,twoN) + 1;
            hy = Int(floor(mod(hind-1,twoN*twoN) / (twoN))) + 1;
            hz = Int(floor((hind-1) / (twoN*twoN))) + 1;
            
            if hx <= griddim[1] && hy <= griddim[2] && hz <= griddim[3]
                gind = hx + griddim[1]*((hy-1) + griddim[2]*(hz-1));
                subrorder[subind] = gind;
                subind += 1;
            end
        end
        return subrorder;
    end
end

# Use this to reorder the nodes in a given grid.
function reorder_grid_recursive!(grid, griddim, type)
    dim = size(grid.allnodes,1);
    
    # get the ordering
    rordering = get_recursive_order(type, dim, griddim);
    # invert it for transfering the grid
    rordering = invert_ordering(rordering);
    
    reorder_grid_nodes!(grid, rordering);
    return grid;
end


#####################################################################################
# ordering rule sets

# Build a Morton Ordering
function morton_ordering_rules(dim)
    if dim == 1
        return RecursiveOrdering(1,[[1]],[[1]]);
    elseif dim == 2
        return RecursiveOrdering(1,[[1,1,1,1]], [[1,2,3,4]]);
    elseif dim == 3
        return RecursiveOrdering(1,[[1,1,1,1,1,1,1,1]], [[1,2,3,4,5,6,7,8]]);
    end
end

# Build a Hilbert Ordering
function hilbert_ordering_rules(dim)
    if dim == 1
        return RecursiveOrdering(1,[[1]], [[1]]);
    elseif dim == 2
        return RecursiveOrdering(4,
                                [[2,3,1,1], [1,2,4,2], [3,1,3,4], [4,4,2,3]], # state indices
                                [[1,4,2,3], [1,2,4,3], [3,4,2,1], [3,2,4,1]]);# indexing order
    elseif dim == 3
        # Oh boy rotations for each state
        rots = [[],
                [1,2],
                [1,2,1,2],
                [2],
                [1,2,2],
                [3],
                [1,1,1],
                [1,1,2],
                [3,3,3],
                [2,3],
                [2,1,1,1],
                [3,3],
                [2,2,2],
                [1],
                [1,1,3],
                [1,2,2,2],
                [1,3],
                [2,2],
                [1,3,3,3],
                [3,2,2,1],
                [1,1],
                [3,3,2],
                [2,2,1],
                [2,2,3]];
        first = [1,4,2,3,8,5,7,6];
        orders = [];
        for i=1:24
            push!(orders, rotations(first, rots[i]));
        end
        
        # state rules
        rules = [[2,12,3,3,20,12,17,17],
                [3,1,11,21,16,1,16,21],
                [1,18,19,19,2,10,2,10],
                [], [],[], [],[],[],  # These ones aren't used, so I'll just leave them empty
                [18,3,12,11,18,20,12,20],
                [17,17,21,12,2,10,2,10],
                [11,11,1,10,19,19,1,16],
                [], [],[],
                [18,2,12,2,18,17,12,19],
                [20,16,20,16,1,18,11,11],
                [21,10,3,3,21,16,17,17],
                [20,16,20,16,3,3,21,12],
                [10,1,10,21,17,1,19,21],
                [11,11,2,18,19,19,20,18],
                [], [],[]];
        
        return RecursiveOrdering(24, rules, orders);
    end
end

function rotation_3d(a,axis)
    r = [5 6 1 2 7 8 3 4 ; 2 6 4 8 1 5 3 7 ; 2 4 1 3 6 8 5 7];
    b = zeros(Int,length(a));
    for i=1:length(a)
        b[r[axis,i]] = a[i];
    end
    return b;
end

function rotations(a,rots)
    b = a;
    for i=1:length(rots)
        b = rotation_3d(b,rots[i]);
    end
    return b;
end
############################################################################################

# Build the ordering for a (2^N)^D grid
function build_ordering_2d(ord, N, invert = true)
    nnodes = (2^N)^2;
    lexorder = zeros(Int,nnodes,2);
    twoN = 2^N;
    for j=1:twoN
        for i=1:twoN
            ni = i + twoN*(j-1);
            lexorder[ni,1] = i;
            lexorder[ni,2] = j;
        end
    end
    
    bbox = zeros(4);
    bbox[1] = 1;
    bbox[3] = 1;
    bbox[2] = twoN;
    bbox[4] = twoN;
    
    center = [(bbox[1]+bbox[2])/2, (bbox[3]+bbox[4])/2];
    tmpbbox = zeros(4);
    tmpcenter = zeros(2);
    
    result = zeros(Int64, nnodes);
    for ni=1:nnodes
        vertex = lexorder[ni,:];
        for bi=1:4
            tmpbbox[bi] = bbox[bi];
        end
        tmpcenter[1] = center[1];
        tmpcenter[2] = center[2];
        
        # Find the index for this node
        index = Int64(0);
        state = 1;
        for i=1:N
            stupid = index;
            stupid = state;
            index = index<<2;
            # Which quadrant is it in?
            quad = 0;
            if vertex[1] > tmpcenter[1]
                quad = quad+1;
            end
            if vertex[2] > tmpcenter[2]
                quad = quad+2;
            end
            quad += 1; # 1-based index
            #println("ni="*string(ni)*" level="*string(i)*" quad="*string(quad));
            
            # Add this quad's ordering position to index 
            index += ord.orders[state][quad] - 1;
            
            # Update the state to that quad's state according to rule for this one.
            state = ord.rules[state][quad];
            
            # Shrink the box to one quad
            for j=1:2
                if vertex[j] > tmpcenter[j]
                    tmpbbox[2*j-1] = tmpcenter[j];
                else
                    tmpbbox[2*j] = tmpcenter[j];
                end
            end
            tmpcenter = [(tmpbbox[1]+tmpbbox[2])/2, (tmpbbox[3]+tmpbbox[4])/2];
        end
        
        # add one for 1-based index
        index += 1;
        #println("index for ni="*string(ni)*" is "*string(index));
        if invert
            result[ni] = index;
        else
            result[index] = ni;
        end
    end
    
    return result;
end

# Build the ordering for a (2^N)^D grid
function build_ordering_3d(ord, N, invert = true)
    nnodes = (2^N)^3;
    lexorder = zeros(Int,nnodes,3);
    twoN = 2^N;
    for k=1:twoN
        for j=1:twoN
            for i=1:twoN
                ni = i + twoN*((j-1) + twoN*(k-1));
                lexorder[ni,1] = i;
                lexorder[ni,2] = j;
                lexorder[ni,3] = k;
            end
        end
    end
    
    bbox = zeros(6);
    bbox[1] = 1;
    bbox[3] = 1;
    bbox[5] = 1;
    bbox[2] = twoN;
    bbox[4] = twoN;
    bbox[6] = twoN;
    
    center = [(bbox[1]+bbox[2])/2, (bbox[3]+bbox[4])/2, (bbox[5]+bbox[6])/2];
    tmpbbox = zeros(6);
    tmpcenter = zeros(3);
    
    result = zeros(Int64, nnodes);
    for ni=1:nnodes
        vertex = lexorder[ni,:];
        for bi=1:6
            tmpbbox[bi] = bbox[bi];
        end
        tmpcenter[1] = center[1];
        tmpcenter[2] = center[2];
        tmpcenter[3] = center[3];
        
        # Find the index for this node
        index = Int64(0);
        state = 1;
        for i=1:N
            stupid = index;
            stupid = state;
            index = index<<3;
            # Which octant is it in?
            octant = 0;
            if vertex[1] > tmpcenter[1]
                octant = octant+1;
            end
            if vertex[2] > tmpcenter[2]
                octant = octant+2;
            end
            if vertex[3] > tmpcenter[3]
                octant = octant+4;
            end
            octant += 1; # 1-based index
            #println("ni="*string(ni)*" level="*string(i)*" oct="*string(octant));
            
            # Add this octant's ordering position to index 
            index += ord.orders[state][octant] - 1;
            
            # Update the state to that octant's state according to rule for this one.
            state = ord.rules[state][octant];
            
            # Shrink the box to one octant
            for j=1:3
                if vertex[j] > tmpcenter[j]
                    tmpbbox[2*j-1] = tmpcenter[j];
                else
                    tmpbbox[2*j] = tmpcenter[j];
                end
            end
            tmpcenter = [(tmpbbox[1]+tmpbbox[2])/2, (tmpbbox[3]+tmpbbox[4])/2, (tmpbbox[5]+tmpbbox[6])/2];
        end
        
        # add one for 1-based index
        index += 1;
        #println("index for ni="*string(ni)*" is "*string(index));
        if invert
            result[ni] = index;
        else
            result[index] = ni;
        end
    end
    
    return result;
end

function invert_ordering(ord)
    iord = zeros(Int,length(ord));
    for i=1:length(ord)
        iord[ord[i]] = i;
    end
    return iord;
end