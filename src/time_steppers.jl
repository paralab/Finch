#=
# A set of time steppers
=#

# mutable struct Stepper
#     type::String;       # The constant for the stepper type
#     implicit::Bool;     # implicit or explicit
#     Nsteps::Int;        # number of steps
#     dt::Float64;        # step size
#     cfl::Float64;       # CFL number
#     stages::Int;        # how many stages
#     a::Array{Float64};  # for RK steppers
#     b::Array{Float64};  # for RK steppers
#     c::Array{Float64};  # for RK steppers
    
#     Stepper(t, c) = new(t, false, 0, 0, c, 0, [], [], []);
# end

function init_stepper(dx, stepper::Stepper)
    if stepper.type == EULER_EXPLICIT
        stepper.implicit = false;
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == EULER_IMPLICIT
        stepper.implicit = true;
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == PECE
        stepper.implicit = false;
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == CRANK_NICHOLSON
        stepper.implicit = true;
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == LSRK4
        stepper.implicit = false;
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 5;
        stepper.a = [0, -0.41789047449985195, -1.192151694642677, -1.6977846924715279, -1.5141834442571558];
        stepper.b = [0.14965902199922912, 0.37921031299962726,  0.8229550293869817, 0.6994504559491221, 0.15305724796815198];
        stepper.c = [0.0, 0.14965902199922912, 0.37040095736420475, 0.6222557631344432, 0.9582821306746903];
        
        return stepper;
        
    elseif stepper.type == RK4
        stepper.implicit = false;
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        stepper.stages = 4;
        stepper.a = [   0   0   0   0;
                        0.5 0   0   0;
                        0   0.5 0   0;
                        0   0   1   0];
        stepper.b = [1/6, 1/3, 1/3, 1/6];
        stepper.c = [0.0, 0.5, 0.5, 1];
        
        return stepper;
        
    elseif stepper.type == "custom multistep"
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        T = finch_state.prob.end_time;
        stepper.dt = stepper.cfl*dx;
        stepper.Nsteps = ceil(T/stepper.dt);
        stepper.dt = T/stepper.Nsteps;
        # stages, a, b, c have been set already
        
        return stepper;
    end
end
