#=
Functions to reformat symbolic layer expressions for time steppers.
This is done in different ways depending on stepper and discretization
=#

function reformat_for_stepper(stepper, lhs, rhs, face_lhs, face_rhs, nlt)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    newfacelhs = [];
    newfacerhs = [];
    newnlt = [];
    if face_lhs === nothing
        has_surface = false;
        face_lhs = [];
        face_rhs = [];
        for vi=1:length(rhs)
            push!(face_lhs, []);
            push!(face_rhs, []);
        end
    else
        has_surface = true;
    end
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        newfacelhs = copy(face_rhs);
        newfacerhs = copy(face_rhs);
        newnlt = copy(nlt);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                if has_surface
                    (newlhs[vi], newrhs[vi], newfacelhs[vi], newfacerhs[vi], newnlt[vi]) = 
                        reformat_for_stepper(stepper, (lhs[1][vi], lhs[2][vi]), rhs[vi], face_lhs[vi], face_rhs[vi], nlt[vi]);
                else
                    (newlhs[vi], newrhs[vi], newnlt[vi]) = 
                        reformat_for_stepper(stepper, (lhs[1][vi], lhs[2][vi]), rhs[vi], nothing, nothing, nlt[vi]);
                end
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
                newfacelhs[vi] = face_lhs[vi];
                newfacerhs[vi] = face_rhs[vi];
                newnlt[vi] = nlt[vi];
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
            for i=1:length(nlt)
                nlt[i] = nlt[i]*dt; # dt*nlt
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs[2]); # lhs1 + dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs + lhs1
            newfacerhs = copy(face_rhs);# dt*facerhs
            newfacelhs = copy(face_lhs);# dt*facelhs
            newnlt = copy(nlt);
            
        elseif stepper == EULER_EXPLICIT || stepper == PECE # lhs1 = lhs1 + dt*rhs - dt*lhs2
            # This version includes the lhs1 + dt*(  ) part
            for i=1:length(lhs[2])
                lhs[2][i] = -lhs[2][i]*dt; # -dt*lhs2
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs)
                face_rhs[i] = face_rhs[i]*dt; # dt*facerhs
            end
            for i=1:length(face_lhs)
                face_lhs[i] = -face_lhs[i]*dt; # -dt*facelhs
            end
            for i=1:length(nlt)
                nlt[i] = nlt[i]*dt; # dt*nlt
            end
            
            newlhs = copy(lhs[1]);# lhs1
            newrhs = copy(rhs);
            append!(newrhs, lhs[2]);# dt*rhs - dt*lhs2
            append!(newrhs, lhs[1]);# dt*rhs - dt*lhs2 + lhs1
            newfacerhs = copy(face_rhs);
            append!(newfacerhs, face_lhs);# dt*facerhs - dt*facelhs
            newfacelhs = [];
            newnlt = copy(nlt);
            
            # # This version does not include the lhs1 + dt*(  ) part
            # for i=1:length(lhs[2])
            #     lhs[2][i] = -lhs[2][i]; # -lhs2
            # end
            # for i=1:length(face_lhs)
            #     face_lhs[i] = -face_lhs[i]; # -facelhs
            # end
            
            # newlhs = copy(lhs[1]);# lhs1
            # newrhs = copy(rhs);
            # append!(newrhs, lhs[2]);# rhs - lhs2
            # newfacerhs = copy(face_rhs);
            # append!(newfacerhs, face_lhs);# facerhs - facelhs
            # newfacelhs = [];
            
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
            for i=1:length(nlt)
                nlt[i] = nlt[i]*dt; # dt*nlt
            end
            
            newlhs = copy(lhs[1]);
            append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            newrhs = copy(rhs);
            append!(newrhs, lhs[1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            append!(newrhs, lhs2r);
            newfacelhs = facelhs2l;
            newfacerhs = face_rhs[1];
            append!(newfacerhs, facelhs2r);
            newnlt = copy(nlt);
            
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
            newnlt = copy(nlt);
            
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
            newnlt = copy(nlt);
            
        end
    end
    
    if has_surface
        return (newlhs, newrhs, newfacelhs, newfacerhs, newnlt);
    else
        return (newlhs, newrhs, newnlt);
    end
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv(stepper, flhs, frhs, slhs, srhs)
    (newflhs, newfrhs) = reformat_for_stepper_fv_flux(stepper, flhs, frhs);
    (newslhs, newsrhs) = reformat_for_stepper_fv_source(stepper, slhs, srhs);
    
    return (newflhs, newfrhs, newslhs, newsrhs);
end

# Special version for FV. Flux term only
function reformat_for_stepper_fv_flux(stepper, lhs, rhs)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs); # copy just to set up the arrays right
        newrhs = copy(rhs);
        for vi=1:length(rhs);
            (newlhs[vi], newrhs[vi]) = reformat_for_stepper_fv_flux(stepper, lhs[vi], rhs[vi]);
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
            # du/dt + lhs = rhs  ->  u(+) + dt*lhs = u(-) + dt*rhs
            # Do the *dt part here, but adding u will be done by the solver
            for i=1:length(lhs)
                lhs[i] = lhs[i]*dt;
            end
            for i=1:length(rhs)
                rhs[i] = rhs[i]*dt;
            end
            
            newlhs = copy(lhs);
            newrhs = copy(rhs);
            
        elseif stepper == CRANK_NICHOLSON 
            printerr("TODO: crank nicholson for FV. Choose something else, please.", fatal=true);
        end
    end
    
    return (newlhs, newrhs);
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv_source(stepper, lhs, rhs)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
        newlhs = copy(rhs); # copy just to set up the arrays right
        newrhs = copy(rhs);
        for vi=1:length(rhs);
            (newlhs[vi], newrhs[vi]) = reformat_for_stepper_fv_source(stepper, lhs[vi], rhs[vi]);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4 || stepper == PECE
            # du/dt = -lhs + rhs
            # for i=1:length(lhs)
            #     lhs[i] = -lhs[i]; # don't change sign here because it is done by parser?
            # end
            for i=1:length(rhs)
                rhs[i] = -rhs[i]; # need to change sign here?
            end
            
            newlhs = []; # put everything in rhs for explicit steppers
            newrhs = copy(rhs);
            append!(newrhs, lhs);
            
        elseif stepper == EULER_IMPLICIT
            # du/dt + lhs = rhs  ->  u(+) + dt*lhs = u(-) + dt*rhs
            # Do the *dt part here, but adding u will be done by the solver
            for i=1:length(lhs)
                lhs[i] = -lhs[i]*dt;
            end
            for i=1:length(rhs)
                rhs[i] = -rhs[i]*dt;
            end
            
            newlhs = copy(lhs);
            newrhs = copy(rhs);
            
        elseif stepper == CRANK_NICHOLSON 
            printerr("TODO: crank nicholson for FV. Choose something else, please.", fatal=true);
        end
    end
    
    return (newlhs, newrhs);
end

# function reformat_for_stepper(stepper, lhs, rhs, face_lhs=nothing, face_rhs=nothing)
#     # rebuild the expressions depending on type of time stepper
#     dt = symbols("dt");
#     newlhs = [];
#     newrhs = [];
    
#     if typeof(rhs) <: Array && length(rhs) > 0 && typeof(rhs[1]) <: Array # recursively work on subarrays
#         newlhs = copy(rhs);
#         newrhs = copy(rhs);
#         for vi=1:length(rhs)
#             if length(lhs[1][vi]) > 0 # this dof has a time derivative term
#                 (newlhs[vi], newrhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], stepper);
#             else # no time derivative for this dof
#                 newlhs[vi] = lhs[2][vi];
#                 newrhs[vi] = rhs[vi];
#             end
#         end
#     else
#         # reformat depending on stepper type
#         if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
#             for i=1:length(lhs[2])
#                 lhs[2][i] = lhs[2][i]*dt; # dt*lhs2
#             end
#             for i=1:length(rhs)
#                 rhs[i] = rhs[i]*dt; # dt*rhs
#             end
            
#             newlhs = copy(lhs[1]);
#             append!(newlhs, lhs[2]); # lhs1 + dt*lhs2
#             newrhs = copy(rhs);
#             append!(newrhs, lhs[1]);# dt*rhs + lhs1
            
#         elseif stepper == EULER_EXPLICIT || stepper == PECE # lhs1 = dt*rhs - dt*lhs2 + lhs1
#             # This version includes the lhs1 + dt*(  ) part
#             for i=1:length(lhs[2])
#                 lhs[2][i] = -lhs[2][i]*dt; # -dt*lhs2
#             end
#             for i=1:length(rhs)
#                 rhs[i] = rhs[i]*dt; # dt*rhs
#             end
            
#             newlhs = copy(lhs[1]);# lhs1
#             newrhs = copy(rhs);
#             append!(newrhs, lhs[2]);# dt*rhs - dt*lhs2
#             append!(newrhs, lhs[1]);# dt*rhs - dt*lhs2 + lhs1
            
#             # # This version doesn't include the lhs1 + dt*(  ) part
#             # for i=1:length(lhs[2])
#             #     lhs[2][i] = -lhs[2][i]; # -lhs2
#             # end
            
#             # newlhs = copy(lhs[1]);# lhs1
#             # newrhs = copy(rhs);
#             # append!(newrhs, lhs[2]);# rhs - lhs2
            
#         elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
#             lhs2l = copy(lhs[2]);
#             lhs2r = copy(lhs[2]);
#             for i=1:length(lhs[2])
#                 lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
#                 lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
#             end
#             for i=1:length(rhs)
#                 rhs[i] = rhs[i]*dt; # dt*rhs
#             end
            
#             newlhs = copy(lhs[1]);
#             append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
#             newrhs = copy(rhs);
#             append!(newrhs, lhs[1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
#             append!(newrhs, lhs2r);
            
#         elseif stepper == LSRK4 # (lhs1) : rhs - lhs2
#             for i=1:length(lhs[2])
#                 lhs[2][i] = -lhs[2][i]; # -lhs2
#             end
#             # for i=1:length(rhs[1])
#             #     rhs[1][i] = rhs[1][i]; # rhs
#             # end
            
#             newlhs = copy(lhs[1]);# lhs1
#             newrhs = copy(rhs);
#             append!(newrhs, lhs[2]);# rhs - lhs2
            
#         elseif stepper == RK4 # (lhs1) : rhs - lhs2
#             for i=1:length(lhs[2])
#                 lhs[2][i] = -lhs[2][i]; # -lhs2
#             end
#             # for i=1:length(rhs[1])
#             #     rhs[1][i] = rhs[1][i]; # rhs
#             # end
            
#             newlhs = copy(lhs[1]);# lhs1
#             newrhs = copy(rhs);
#             append!(newrhs, lhs[2]);# rhs - lhs2
            
#         end
#     end
    
#     return (newlhs, newrhs);
# end