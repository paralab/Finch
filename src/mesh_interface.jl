## Mesh interface functions that we think, which should be supported
# to plugin Mesh library to finch
# this is mostly influenced by the dendro loop strucutres. But We can modify this as we go along. 

# dimention of the simplex
@enum SIMPLEX_DIM begin
    D0=0,
    D1=1,
    D2=2,
    D3=3,
    D4=4
end

# element (support upto 4D elements) loop types
@enum SIMPLEX4_LOOP_TYPE begin
    ALL         =0,
    LOCAL       =1,
    HALO        =2,
    INDEP       =3,
    DEP         =4,
    WRITABLE    =5 
end 

# simplex of dim 3 
@enum SIMPLEX3_LOOP_TYPE begin
    ALL         =0,
    LOCAL       =1,
    HALO        =2,
    INDEP       =3,
    DEP         =4,
    WRITABLE    =5,
    BOUNDRY     =6
end

# simplex of dim 2 
@enum SIMPLEX2_LOOP_TYPE begin
    ALL         =0,
    LOCAL       =1,
    HALO        =2,
    INDEP       =3,
    DEP         =4,
    WRITABLE    =5, 
    BOUNDRY     =6
end

# simplex of dim 1 
@enum SIMPLEX1_LOOP_TYPE begin
    ALL         =0,
    LOCAL       =1,
    HALO        =2,
    INDEP       =3,
    DEP         =4,
    WRITABLE    =5, 
    BOUNDRY     =6
end

# simplex of dim 0 (vertices) 
@enum SIMPLEX0_LOOP_TYPE begin
    ALL         =0,
    LOCAL       =1,
    HALO        =2,
    INDEP       =3,
    DEP         =4,
    WRITABLE    =5, 
    BOUNDRY     =6
end


# init function to initialize the iterator loop. 
function _init{LoopType}() where {LoopType}
end

# increment to the next iterator. 
function _next{LoopType}() where {LoopType}
end

# returns the end iterator
function _end{LoopType}() where {LoopType}
end

# get the current iterator. 
function _curr{LoopType}() where {LoopType}
end

# get the indices of the specifid iterator. 
# curr: current simplex iterator. 
# dim: get the indices of simplex dim for current simplex
# example: if you are iterating over D3 simplex, 
# you can get the indices of D2(surface), D1(edges), D0(vertices)

#= 
function get_indices{Iter}(Iter curr, SIMPLEX_DIM dim) where {Iter}
#end

# dim : dim of the geometric object corresponding to the vector (4d, 3d, 2d, 1d, 0d),
       #if it is 0 then it will be nodal vector, 1 then edge vector, 2 surface vector etc. 
# with_halo : allocate mem. for halo regions as well. 
# dof : degrees of freedoms
function alloc_vector{T}(SIMPLEX_DIM dim, with_halo=false, dof=1) where {T}
end

# deallocation for the vector. 
function dealloc_vector{Vec}(Vec vec) where {Vec}
end

# perform the halo sync. between proc. for a given vector (async begin non blocking call).
function sync_vector_begin{Vec}(Vec vec, with_halo=false, dof=1) where {Vec}
end 

# perform the halo sync. between proc. for a given vector (async end non blocking call).
function sync_vector_end{Vec}(Vec vec, with_halo=false, dof=1) where {Vec}
end

=#