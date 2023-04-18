import os
path = os.path.dirname(os.path.realpath(__file__))
folders = []
mpi = 0
for i in sorted(os.listdir(path)):
  if os.path.isdir(os.path.join(path, i)) and 'results_' in i and 'dendro' not in i:
    folders.append(i)


f = open("timestep.pvd", "w");
f.write("<?xml version=\"1.0\"?>\n")
f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
f.write("<Collection>\n");

filename = "ns" #may need to change this


for i in range (0,len(folders)):
    timestep = filename + folders[i][-6:] + ".pvtu"
    fname = folders[i] + "/" + timestep;
    f.write("<DataSet timestep=\"" + str(i) + "\" group=\"\" part=\"0\" file=\"" + fname + "\"/>\n")

f.write("</Collection>\n")
f.write("</VTKFile>\n")
f.close()
