#=
Macros for the interface.
=#
export @preStepFunction, @postStepFunction, @callbackFunction

#= Assume something like this
myf = @callbackFunction(
    function f(a,b,c)
        body...
    end
)
=#
macro callbackFunction(f)
    fname = string(f.args[1].args[1]);
    fargs = [string(f.args[1].args[i]) for i=2:length(f.args[1].args)];
    fbody = string(f.args[2]);
    return esc(:(callbackFunction($f, name=$fname, args=$fargs, body=$fbody)));
end

macro preStepFunction(f)
    return esc(quote
        function PRE_STEP_FUNCTION()
            $f
        end
        preStepFunction(PRE_STEP_FUNCTION);
    end)
end

macro postStepFunction(f)
    return esc(quote
        function POST_STEP_FUNCTION()
            $f
        end
        postStepFunction(POST_STEP_FUNCTION);
    end)
end
