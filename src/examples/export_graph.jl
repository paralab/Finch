### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

#=
Import or create a mesh, then export the graph info.
That is all this script does.

An edge exists between two elements that share a node.
The weight is how many shared nodes.
The format is 
- number of dimensions (2 or 3)
- number of elements (integer)
- number of edges (integer)
- element center coordinates (two or three floats each)
- edge pairs and weight (three integers each)
=#
initFinch("export_graph");


domain(2) # Change the dimensionality. Can be 2 or 3.
mesh("src/examples/utriangle.msh") # Give the .msh or .mesh file (Gmsh or Meddit formats)
exportMeshGraph("graphinfo.txt") # Choose your output file name


finalizeFinch()
