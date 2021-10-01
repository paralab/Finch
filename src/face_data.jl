#=
# A struct containing face info
=#
export FaceData

mutable struct FaceData
    id;                  # face id
    normal;   # normal of this side of this face
    nodes;      	 # nodes on this side of this face
    isbdy;
    FaceData(id, face2v, normals) = (
		(normal, nodes, isbdy) = init(id, face2v, normals);
		new(id,normal,nodes,isbdy);
	)
    FaceData(id,normal,nodes,isbdy) = (
		new(id,normal,nodes,isbdy);
	)
	
    
end

function init(id, face2v, nomals)
	
	nodes = face2v[:,:,id];
	normal = [normals[:,id] -normals[:,id]]; 
	isbdy = 1;
	return [normal, nodes, isbdy];

end

