

n = 6;
ord = 5;

init_finch("reordering");

@domain(3)
@functionSpace(LEGENDRE, ord)

@mesh(HEXMESH, n)

@variable(u)
@variable(q1)
@variable(q2)
@variable(q3)
@variable(q4)
@variable(q5)
@variable(q6)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)
@boundary(q1, 1, DIRICHLET, 0)
@boundary(q2, 1, DIRICHLET, 0)
@boundary(q3, 1, DIRICHLET, 0)
@boundary(q4, 1, DIRICHLET, 0)
@boundary(q5, 1, DIRICHLET, 0)
@boundary(q6, 1, DIRICHLET, 0)

@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

@coefficient(g, "(-3*pi*pi + 5)*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm([q1,q2,q3,q4,q5,q6], ["-dot(grad(q1), grad(v)) + q2*v + q3*v + q4*v + q5*v + q6*v - g*v",
                                "-dot(grad(q2), grad(v)) + q3*v + q4*v + q5*v + q6*v + q1*v - g*v",
                                "-dot(grad(q3), grad(v)) + q4*v + q5*v + q6*v + q1*v + q2*v - g*v",
                                "-dot(grad(q4), grad(v)) + q5*v + q6*v + q1*v + q2*v + q3*v - g*v",
                                "-dot(grad(q5), grad(v)) + q6*v + q1*v + q2*v + q3*v + q4*v - g*v",
                                "-dot(grad(q6), grad(v)) + q1*v + q2*v + q3*v + q4*v + q5*v - g*v"])

times = 2;
timings = zeros(3);