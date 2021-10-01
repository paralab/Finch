# Some utils for DG issues
export dg_nodes_element_to_face, dg_nodes_face_to_element

# For face integrals, change a vector from full element to face components only
# and the other way around just inserting zeros.
function dg_nodes_element_to_face(v, refel, face)
    nfp = refel.Nfp[face];
    vf = zeros(nfp);
    
    for i=1:nfp
        vf[i] = v[refel.face2local[face][i]];
    end
    
    return vf;
end

function dg_nodes_face_to_element(vf, refel, face)
    nfp = refel.Nfp[face];
    v = zeros(refel.Np);
    
    for i=1:nfp
        v[refel.face2local[face][i]] = vf[i];
    end
    
    return v;
end