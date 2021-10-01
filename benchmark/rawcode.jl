# Find a way to run these example in isolation
using Finch


init_finch("poisson1d")

useLog("poisson1dlog")

domain(1)

functionSpace(order=3)
mesh(LINEMESH, elsperdim=20)

u = variable("u")
testSymbol("v")
boundary(u, 1, DIRICHLET, 0)

coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) -
             pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
weakForm(u, "-grad(u)*grad(v) - f*v")

import Finch: grid_data, mesh_data, refel 
dofs_per_node = 1

N1 = size(grid_data.allnodes,2);
Nn = dofs_per_node * N1;
Np = refel.Np;
nel = mesh_data.nel;

rhsvec = zeros(Nn);
lhsmatI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
lhsmatJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
lhsmatV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
allocated_vecs = [rhsvec, lhsmatI, lhsmatJ, lhsmatV];

Finch.CGSolver.assemble(u, Finch.bilinears[u.index], Finch.linears[u.index], allocated_vecs)

print("END")

#SUITE["poisson1d-assemble"] = @benchmarkable Finch.CGSolver.assemble($u,
                         #$Finch.bilinears[u.index], $Finch.linears[u.index])