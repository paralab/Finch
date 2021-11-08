#=
Represents cells in the neighborhood of a face.
Has:
- cell indices
- cell center coords
- specified variable values (optional)
=#
struct Neighborhood
    cells::Array{Array{Int,1},1}        # indices for cells on left and right: [[L1,L2,...], [R1,R2,...]]
    centers::Array{Array{Float64,2},1}  # cell center coords [[L1x,L2x; L1y,L2y; L1z,L2z], [similar]]
    values::Union{Array{Array{Float64,2},1}, Array{Array{Float64,1},1}, Nothing}   # variable values in those cells or nothing
    
    Neighborhood(ind, center) = new(ind, center, nothing);
    Neighborhood(ind, center, var::Variable) = new(ind, center, get_neighborhood_variable_values(ind, var));
    Neighborhood(ind, center, var::Array{Array{Float64,1},1}) = new(ind, center, var);
    Neighborhood(ind, center, var::Array{Array{Float64,2},1}) = new(ind, center, var);
end

function get_neighborhood_variable_values(cells::Array{Array{Int,1},1}, var::Array{Variable,1})
    nvars = 0;
    for vi=1:length(var)
        nvars += size(var[vi].values,1);
    end
    vals = [zeros(nvars, length(cells[1])), zeros(nvars, length(cells[2]))];
    
    st = 1;
    for vi=1:length(var)
        en = st + size(var[vi].values,1) - 1;
        vals[1][st:en,:] = var[vi].values[:,cells[1]];
        vals[2][st:en,:] = var[vi].values[:,cells[2]];
        st = en+1;
    end
    
    return vals;
end
function get_neighborhood_variable_values(cells::Array{Array{Int,1},1}, var::Variable)
    return [var.values[:,cells[1]], var.values[:,cells[2]]];
end

Base.length(n::Neighborhood) = 1;
