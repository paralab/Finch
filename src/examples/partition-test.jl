
package_installed = false;
if package_installed
    using Finch
else
    if !@isdefined(Finch)
        include("../Finch.jl");
        using .Finch
    end
end

init_finch("partitionTest");
dim = 3;

if dim == 1
    domain(1)
    mesh(LINEMESH, elsperdim=40)
    np = 5;
elseif dim == 2
    domain(2)
    mesh(QUADMESH, elsperdim=40)
    np = 7;
else
    domain(3)
    mesh(HEXMESH, elsperdim=20)
    np = 12;
end

m = Finch.mesh_data;

epart = Finch.get_element_partitions(m, np);

(r0,g0) = Finch.partitioned_grid_from_mesh(m, epart);

rs = [r0];
gs = [g0];

for i=2:np
    Finch.config.proc_rank = i-1;
    (tmpr,tmpg) = Finch.partitioned_grid_from_mesh(m, epart);
    push!(rs, tmpr);
    push!(gs, tmpg);
end

########## Test the partitioned grids for correctness ##############
println("Testing partitions.");

# test for correct number of elements in each
e_count = zeros(Int, np);
for ei=1:length(epart)
    e_count[epart[ei]+1] += 1;
end
println("Expected owned elements: "*string(e_count));

nerror = 0;
for i=1:np
    nel = length(gs[i].element_owner);
    mynel = 0;
    for ei=1:nel
        if gs[i].element_owner[ei] < 0
            mynel += 1;
        end
    end
    println("partition "*string(i-1)*" owns "*string(mynel)*" ("*string(nel)*" including ghosts)")
    if !(e_count[i] == mynel)
        global nerror += 1;
    end
end
if nerror > 0
    println("There were "*string(nerror)*" errors!")
else
    println("All element counts are correct")
end
println("");

# test that ghost owners are correct
nerror = 0;
for i=1:np
    nel = length(gs[i].element_owner);
    for ei=1:nel
        glb = gs[i].grid2mesh[ei];
        owner = gs[i].element_owner[ei];
        if owner < 0
            owner = i-1;
        end
        if !(epart[glb] == owner)
            global nerror += 1;
            println("error: partition "*string(np-1)*" thinks element "*string(glb)*" belongs to "*string(owner)*" but it belongs to "*string(epart[glb]));
        end
    end
end
if nerror == 0
    println("All ghost owners are correct")
end
println("");

# test boundary face counts
bdry_count = zeros(Int, np);
bdry_face = zeros(Int, np);
for i=1:np
    nface = size(gs[i].face2element,2);
    for fi=1:nface
        e1 = gs[i].face2element[1,fi];
        e2 = gs[i].face2element[2,fi];
        if e2 == 0
            bdry_count[i] += 1;
        end
    end
    bdry_face[i] = length(gs[i].bdryface[1]);
end
total_bdry_face = 0;
for fi=1:length(m.bdryID)
    if m.bdryID[fi] > 0
        global total_bdry_face += 1;
    end
end
println(string(total_bdry_face)*" bdry faces in mesh,");
println(string(sum(bdry_face))*" bdry faces in grids("*string(bdry_face)*"),");
println(string(sum(bdry_count))*" faces with an empty neighbor in grids("*string(bdry_count)*")");
if total_bdry_face == sum(bdry_face) && total_bdry_face == sum(bdry_count)
    println("Numbers of boundary faces are correct")
else
    println("Error in numbers of boundary faces")
end

