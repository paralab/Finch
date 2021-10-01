using Finch
using Test

@testset "Finch.jl" begin
    init_finch("poisson1d");

    useLog("poisson1dlog")

    domain(1)                      # dimension
    functionSpace(order=3)         # basis function polynomial order

    mesh(LINEMESH, elsperdim=20)   # build uniform LINEMESH with 20 elements

    u = variable("u")              # make a scalar variable with symbol u
    testSymbol("v")                # sets the symbol for a test function

    boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

    coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
    weakForm(u, "-grad(u)*grad(v) - f*v")

    solve(u);
    
    allerr = zeros(size(Finch.grid_data.allnodes,2));
    exact(x) = sin(10*pi*x)*sin(pi*x);

    for i=1:size(Finch.grid_data.allnodes,2)
        x = Finch.grid_data.allnodes[1,i];
        err = abs(u.values[i] - exact(x));
        allerr[i] = err;
    end
    maxerr = maximum(abs, allerr);
    println("max error = "*string(maxerr));
    
    @test(maxerr < 0.01)
    
    # Write your tests here.
end
