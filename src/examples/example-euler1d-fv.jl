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
n = 50 # number of elements
mesh(LINEMESH, elsperdim=n)

# Set the order for flux
finiteVolumeOrder(2);

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

# Transform between variable types
function U2Q!(U,Q)
    local gamma = 1.4;
    if typeof(U[1]) == Finch.Variable
        Q[1].values .= U[1].values;
        Q[2].values .= U[1].values .* U[2].values;
        Q[3].values .= U[3].values ./ (gamma-1) + 0.5 .* U[1].values .* U[2].values .^2;
    else
        Q[1] = U[1];
        Q[2] = U[1] * U[2];
        Q[3] = U[3] / (gamma-1) + 0.5 * U[1] * U[2] * U[2];
    end
end
function Q2U!(U,Q)
    local gamma = 1.4;
    if typeof(U[1]) == Finch.Variable
        U[1].values .= Q[1].values;
        U[2].values .= Q[2].values ./ Q[1].values;
        U[3].values .= (gamma-1) .* (Q[3].values .- 0.5 .* Q[2].values .* Q[2].values ./ Q[1].values);
    else
        U[1] = Q[1];
        U[2] = Q[2] / Q[1];
        U[3] = (gamma-1) * (Q[3] - 0.5 * Q[2] * Q[2] / Q[1]);
    end
end

# Shock tube initial conditions
initial(r, "x<0.3 ? 1 : 0.125")
initial(u, "x<0.3 ? 0.75 : 0")
initial(p, "x<0.3 ? 1 : 0.1")
Finch.eval_initial_conditions(); # set U initial conditions
U2Q!(U,Q); # set Q initial conditions

# Boundary conditions
boundary(q1, 1, NO_BC)
boundary(q2, 1, NO_BC)
boundary(q3, 1, NO_BC)

gamma = coefficient("gamma", 1.4)

function the_flux_function(U)
    local r = U[1];
    local u = U[2];
    local p = U[3];
    local gamma = 1.4;
    
    FU = [ r * u, 
           p + r * u*u, 
           u * (p + p/(gamma-1) + r*u*u*0.5)];
    
    # This would be generated by some flux splitting operator
    c = sqrt(gamma*p/r)
    M = u/c
    if M < -1
        Fm = FU;
        Fp = zeros(3);
    elseif M > 1
        Fm = zeros(3);
        Fp = FU;
    elseif M >= -1 && M < 0
        Fp = r*(u+c)/(2*gamma) .* [1, (u+c), (u+c)^2/2+(3-gamma)/(gamma-1)*c^2/2];
        Fm = FU - Fp;
    else # 0 < M < 1
        Fm = r*(u-c)/(2*gamma) .* [1, (u-c), (u-c)^2/2+(3-gamma)/(gamma-1)*c^2/2];
        Fp = FU - Fm;
    end
    
    return (Fm, Fp);
end

callbackFunction(the_flux_function);

#@exportCode("fveuler1dcodeout") # uncomment to export generated code to a file
@importCode("fveuler1dcode") # uncomment to import code from a file

function pre_step_fun()
    Q2U!(U,Q);
end
function post_step_fun()
    Q2U!(U,Q);
end
preStepFunction(pre_step_fun);
postStepFunction(post_step_fun);

solve(Q)

finalize_finch()

##### Uncomment below to compare to plot

x = Finch.fv_info.cellCenters[:];
# n = length(x);
# ic = zeros(n);
# for i=1:n
#     # ic[i] = (x[i]>0.1 && x[i]<0.5) ? 1 : 0
#     ic[i] = x[i]<0.5 ? -1 : 1
#     # ic[i] = sin(2*pi*x[i])
# end

using Plots
pyplot();
display(plot([x x x], [r.values[:] u.values[:] p.values[:]], markershape=:circle, label=["r" "u" "p"]))
