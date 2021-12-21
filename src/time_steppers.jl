#=
# A set of time steppers
=#

mutable struct Stepper
    type::String;       # The constant for the stepper type
    implicit::Bool;     # implicit or explicit
    Nsteps::Int;        # number of steps
    dt::Float64;        # step size
    cfl::Float64;       # CFL number
    stages::Int;        # how many stages
    a::Array{Float64};  # for RK steppers
    b::Array{Float64};  # for RK steppers
    c::Array{Float64};  # for RK steppers
    
    Stepper(t, c) = new(t, false, 0, 0, c, 0, [], [], []);
end

function init_stepper(x, stepper)
    if stepper.type == EULER_EXPLICIT
        stepper.implicit = false;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == EULER_IMPLICIT
        stepper.implicit = true;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = (prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == PECE
        stepper.implicit = false;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == CRANK_NICHOLSON
        stepper.implicit = true;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == LSRK4
        stepper.implicit = false;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 5;
        stepper.a = [0; -0.41789047449985195; -1.192151694642677; -1.6977846924715279; -1.5141834442571558];
        stepper.b = [0.14965902199922912; 0.37921031299962726;  0.8229550293869817; 0.6994504559491221; 0.15305724796815198];
        stepper.c = [0.0; 0.14965902199922912; 0.37040095736420475; 0.6222557631344432; 0.9582821306746903];
        
        return stepper;
        
    elseif stepper.type == RK4
        stepper.implicit = false;
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 4;
        stepper.a = [   0   0   0   0;
                        0.5 0   0   0;
                        0   0.5 0   0;
                        0   0   1   0];
        stepper.b = [1/6; 1/3; 1/3; 1/6];
        stepper.c = [0.0; 0.5; 0.5; 1];
        
        return stepper;
    end
end

function reformat_for_stepper(lhs, rhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], stepper);
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
            end
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
            for i=1:length(lhs[2])
                lhs[2][i] = lhs[2][i]*dt; # dt*lhs2
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt; # dt*rhs
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs[2]); # lhs1 + dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs + lhs1
            
        elseif stepper == EULER_EXPLICIT || stepper == PECE # lhs1 = dt*rhs - dt*lhs2 + lhs1
            # # This version includes the lhs1 + dt*(  ) part
            # for i=1:length(lhs[2])
            #     lhs[2][i] = -lhs[2][i]*dt; # -dt*lhs2
            # end
            # for i=1:length(rhs)
            #     rhs[i] = rhs[i]*dt; # dt*rhs
            # end
            
            # newlhs = copy(lhs[1]);# lhs1
            # newrhs = copy(rhs);
            # append!(newrhs, lhs[2]);# dt*rhs - dt*lhs2
            # append!(newrhs, lhs[1]);# dt*rhs - dt*lhs2 + lhs1
            
            # This version doesn't include the lhs1 + dt*(  ) part
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            
        elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
            lhs2l = copy(lhs[2]);
            lhs2r = copy(lhs[2]);
            for i=1:length(lhs[2])
                lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
                lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt; # dt*rhs
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            append!(newrhs, lhs2r);
            
        elseif stepper == LSRK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]; # rhs
            # end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            
        elseif stepper == RK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]; # rhs
            # end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            
        end
    end
    
    return (newlhs, newrhs);
end

function reformat_for_stepper(lhs, rhs, face_lhs, face_rhs,stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    newfacelhs = [];
    newfacerhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        newfacelhs = copy(face_rhs);
        newfacerhs = copy(face_rhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi], newfacelhs[vi], newfacerhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], face_lhs[vi], face_rhs[vi], stepper);
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
                newfacelhs[vi] = face_lhs[vi];
                newfacerhs[vi] = face_rhs[vi];
            end
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
            for i=1:length(lhs[2])
                lhs[2][i] = lhs[2][i]*dt; # dt*lhs2
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs)
                face_rhs[i] = face_rhs[i]*dt; # dt*facelhs
            end
            for i=1:length(face_lhs)
                face_lhs[i] = face_lhs[i]*dt; # dt*facerhs
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs[2]); # lhs1 + dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs + lhs1
            newfacerhs = copy(face_rhs);# dt*facerhs
            newfacelhs = copy(face_lhs);# dt*facelhs
            
        elseif stepper == EULER_EXPLICIT || stepper == PECE # lhs1 = lhs1 + dt*rhs - dt*lhs2
            # # This version includes the lhs1 + dt*(  ) part
            # for i=1:length(lhs[2])
            #     lhs[2][i] = -lhs[2][i]*dt; # -dt*lhs2
            # end
            # for i=1:length(rhs)
            #     rhs[i] = rhs[i]*dt; # dt*rhs
            # end
            # for i=1:length(face_rhs)
            #     face_rhs[i] = face_rhs[i]*dt; # dt*facerhs
            # end
            # for i=1:length(face_lhs)
            #     face_lhs[i] = -face_lhs[i]*dt; # -dt*facelhs
            # end
            
            # newlhs = copy(lhs[1]);# lhs1
            # newrhs = copy(rhs);
            # append!(newrhs, lhs[2]);# dt*rhs - dt*lhs2 + lhs1
            # #append!(newrhs, lhs[1]);
            # newfacerhs = copy(face_rhs);
            # append!(newfacerhs, face_lhs);# dt*facerhs - dt*facelhs
            # newfacelhs = [];
            
            # This version does not include the lhs1 + dt*(  ) part
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            for i=1:length(face_lhs)
                face_lhs[i] = -face_lhs[i]; # -facelhs
            end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            newfacerhs = copy(face_rhs);
            append!(newfacerhs, face_lhs);# facerhs - facelhs
            newfacelhs = [];
            
        elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
            lhs2l = copy(lhs[2]);
            lhs2r = copy(lhs[2]);
            facelhs2l = copy(face_lhs);
            facelhs2r = copy(face_lhs);
            for i=1:length(lhs[2])
                lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
                lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(face_lhs)
                facelhs2l[i] = Basic(0.5)*facelhs2l[i]*dt; # 0.5*dt*lhs2
                facelhs2r[i] = Basic(-0.5)*facelhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs)
                face_rhs[i] = face_rhs[i]*dt; # dt*facelhs
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            append!(newrhs, lhs2r);
            newfacelhs = facelhs2l;
            newfacerhs = face_rhs[1];
            append!(newfacerhs, facelhs2r);
            
        elseif stepper == LSRK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            for i=1:length(face_lhs)
                face_lhs[i] = -face_lhs[i]; # -facelhs
            end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            newfacerhs = copy(face_rhs);
            append!(newfacerhs, face_lhs);# facerhs - facelhs
            newfacelhs = [];
            
        elseif stepper == RK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]; # -lhs2
            end
            for i=1:length(face_lhs)
                face_lhs[i] = -face_lhs[i]; # -facelhs
            end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# rhs - lhs2
            newfacerhs = copy(face_rhs);
            append!(newfacerhs, face_lhs);# facerhs - facelhs
            newfacelhs = [];
            
        end
    end
    
    return (newlhs, newrhs, newfacelhs, newfacerhs);
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv(flhs, frhs, slhs, srhs, stepper)
    (newflhs, newfrhs) = reformat_for_stepper_fv_flux(flhs, frhs, stepper);
    (newslhs, newsrhs) = reformat_for_stepper_fv_source(slhs, srhs, stepper);
    
    return (newflhs, newfrhs, newslhs, newsrhs);
end

# Special version for FV. Flux term only
function reformat_for_stepper_fv_flux(lhs, rhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs); # copy just to set up the arrays right
        newrhs = copy(rhs);
        for vi=1:length(rhs);
            (newlhs[vi], newrhs[vi]) = reformat_for_stepper_fv_flux(lhs[vi], rhs[vi], stepper);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4 || stepper == PECE
            # du/dt = -lhs + rhs
            for i=1:length(lhs)
                lhs[i] = -lhs[i];
            end
            # for i=1:length(rhs)
            #     rhs[i] = -rhs[i]; # don't change sign for rhs, because already done by parser
            # end
            
            newlhs = []; # put everything in rhs for explicit steppers
            newrhs = copy(rhs);
            append!(newrhs, lhs);
            
        elseif stepper == EULER_IMPLICIT
            
        elseif stepper == CRANK_NICHOLSON 
            
        end
    end
    
    return (newlhs, newrhs);
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv_source(lhs, rhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs); # copy just to set up the arrays right
        newrhs = copy(rhs);
        for vi=1:length(rhs);
            (newlhs[vi], newrhs[vi]) = reformat_for_stepper_fv_source(lhs[vi], rhs[vi], stepper);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4 || stepper == PECE
            # du/dt = -lhs + rhs
            for i=1:length(lhs)
                lhs[i] = -lhs[i];
            end
            for i=1:length(rhs)
                rhs[i] = -rhs[i]; # need to change sign here?
            end
            
            newlhs = []; # put everything in rhs for explicit steppers
            newrhs = copy(rhs);
            append!(newrhs, lhs);
            
        elseif stepper == EULER_IMPLICIT
            
        elseif stepper == CRANK_NICHOLSON 
            
        end
    end
    
    return (newlhs, newrhs);
end