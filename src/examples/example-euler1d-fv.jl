if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("FVeuler1d");

useLog("FVeuler1dlog", level=3)

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4)
timeInterval(0.2)

# Mesh
n = 150 # number of elements
mesh(LINEMESH, elsperdim=n)

# Primitive variables
r = variable("r", SCALAR, CELL)
u = variable("u", SCALAR, CELL)
p = variable("p", SCALAR, CELL)
U = [r, u, p];

# Conserved variables
q1 = variable("q1", SCALAR, CELL)
q2 = variable("q2", SCALAR, CELL)
q3 = variable("q3", SCALAR, CELL)
Q = [q1, q2, q3];

U2Q = variableTransform(U, Q, 
(U) -> (
    gamma = 1.4;
    Q = [U[1],
         U[1] * U[2],
         U[3] / (gamma-1) + 0.5 * U[1] * U[2] * U[2]];
    return Q;
))
Q2U = variableTransform(Q, U, 
(Q) -> (
    gamma = 1.4;
    U = [Q[1],
         Q[2] / Q[1],
         (gamma-1) * (Q[3] - 0.5 * Q[2] * Q[2] / Q[1])];
    return U;
))

# Shock tube initial conditions
initial(r, "x<0.3 ? 1 : 0.125")
initial(u, "x<0.3 ? 0.75 : 0")
initial(p, "x<0.3 ? 1 : 0.1")
Finch.eval_initial_conditions(); # set U initial conditions
transformVariable(U2Q); # set Q initial conditions

# Boundary conditions
boundary(q1, 1, NO_BC)
boundary(q2, 1, NO_BC)
boundary(q3, 1, NO_BC)

gamma = coefficient("gamma", 1.4)

function flux_u2(r1,u1,p1, r2,u2,p2, comp, normal)
    local gamma = 1.4;
    # The normal points from cell 1 to cell 2
    # So "left" is currently cell 1
    if normal < 0
        # swap left and right
        tmpr = r1; tmpu = u1; tmpp = p1;
        r1 = r2;   u1 = u2;   p1 = p2;
        r2 = tmpr; u2 = tmpu; p2 = tmpp;
    end
    c1 = sqrt(gamma*p1/r1);
    M1 = u1/c1;
    c2 = sqrt(gamma*p2/r2);
    M2 = u2/c2;
    
    if comp == 1
        if M1 < -1
            Fp = 0;
        elseif M1 > 1
            Fp = r1 * u1;
        elseif M1 >= -1 && M1 < 0
            Fp = r1*(u1+c1)/(2*gamma);
        else # 0 < M < 1
            Fp = r1*u1 - r1*(u1-c1)/(2*gamma);
        end
        
        if M2 < -1
            Fm = r2*u2;
        elseif M2 > 1
            Fm = 0;
        elseif M2 >= -1 && M2 < 0
            Fm = r2*u2 - r2*(u2+c2)/(2*gamma);
        else # 0 < M < 1
            Fm = r2*(u2-c2)/(2*gamma);
        end
        
    elseif comp == 2
        if M1 < -1
            Fp = 0;
        elseif M1 > 1
            Fp = p1 + r1 * u1*u1;
        elseif M1 >= -1 && M1 < 0
            Fp = r1*(u1+c1)/(2*gamma) * (u1+c1);
        else # 0 < M < 1
            Fp = p1 + r1 * u1*u1 - r1*(u1-c1)/(2*gamma) * (u1-c1);
        end
        
        if M2 < -1
            Fm = p2 + r2 * u2*u2;
        elseif M2 > 1
            Fm = 0;
        elseif M2 >= -1 && M2 < 0
            Fm = p2 + r2 * u2*u2 - r2*(u2+c2)/(2*gamma) * (u2+c2);
        else # 0 < M < 1
            Fm = r2*(u2-c2)/(2*gamma) * (u2-c2);
        end
        
    elseif comp == 3
        if M1 < -1
            Fp = 0;
        elseif M1 > 1
            Fp = u1 * (p1 + p1/(gamma-1) + r1*u1*u1*0.5);
        elseif M1 >= -1 && M1 < 0
            Fp = r1*(u1+c1)/(2*gamma) * ((u1+c1)^2/2+(3-gamma)/(gamma-1)*c1^2/2);
        else # 0 < M < 1
            Fp = u1 * (p1 + p1/(gamma-1) + r1*u1*u1*0.5) - (r1*(u1-c1)/(2*gamma) * ((u1-c1)^2/2+(3-gamma)/(gamma-1)*c1^2/2));
        end
        
        if M2 < -1
            Fm = u2 * (p2 + p2/(gamma-1) + r2*u2*u2*0.5);
        elseif M2 > 1
            Fm = 0;
        elseif M2 >= -1 && M2 < 0
            Fm = u2 * (p2 + p2/(gamma-1) + r2*u2*u2*0.5) - (r2*(u2+c2)/(2*gamma) * ((u2+c2)^2/2+(3-gamma)/(gamma-1)*c2^2/2));
        else # 0 < M < 1
            Fm = r2*(u2-c2)/(2*gamma) * ((u2-c2)^2/2+(3-gamma)/(gamma-1)*c2^2/2);
        end
    end
    return (Fp + Fm) * normal;
end

callbackFunction(flux_u2);

flux(Q, ["flux_u2(left(r),left(u),left(p), right(r),right(u),right(p), $i, normal())" for i in 1:3])

exportCode("fveuler1dcodeout")

@postStepFunction(
    transformVariable(Q2U)
);

solve(Q)

finalize_finch()

##### Uncomment below to compare to plot

x = Finch.fv_info.cellCenters[:];

using Plots
pyplot();
display(plot([x x x], [r.values[:] u.values[:] p.values[:]], markershape=:circle, label=["r" "u" "p"]))
