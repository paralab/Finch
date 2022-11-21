#=
# Optional log writing
=#
export log_entry, log_dump_config, log_dump_prob

# verbosity levels
# 0 = minimal, not very useful
# 1 = basic messages about steps to show progress
# 2 = default, details about each step
# 3 = debug, lots of info for debugging that is hard to read
function set_log_level(state::FinchState, lev)
    state.log_level = lev;
end

# Prints an error to the log as well as the output and if fatal, exits with code 1
function printerr(state::FinchState, msg; fatal=false)
    log_entry(state, "Error: "*msg, 0);
    if state.config.proc_rank == 0
        println("Error: "*msg);
    end
    if fatal
        println("Fatal error on proc "*string(state.config.proc_rank)*", exiting")
        exit(1);
    end
end
function printerr(msg; fatal=false)
    printerr(finch_state, msg, fatal=fatal);
end

function init_log(state::FinchState, name, dir, level=2)
    state.log_file = dir*"/"*name*".txt";
    state.use_log = true;
    state.log_level = level;
    if state.config.proc_rank == 0
        file = open(state.log_file, "w");
        println(file, "######################################");
        println(file, "# Finch Log for: "*state.project_name);
        println(file, "######################################");
        println(file, "(verbosity = "*string(state.log_level)*")");
        close(file)
    end
end

function log_entry(state::FinchState, text, level=2)
    if state.config.proc_rank == 0 && state.use_log && (level <= state.log_level)
        file = open(state.log_file, "a");
        println(file, string(state.log_line_index)*".\t"*text);
        state.log_line_index += 1;
        close(file);
    end
end
function log_entry(text, level=2)
    log_entry(finch_state, text, level);
end

function log_dump_config(state::FinchState)
    if state.use_log && state.config.proc_rank == 0
        file = open(state.log_file, "a");
        println(file, string(state.log_line_index)*".\tDumping configuration:");
        state.log_line_index += 1;
        for f in fieldnames(FinchConfig)
            println(file, string(state.log_line_index)*".\t\t"*string(f)*" = "*string(getfield(c, f)));
            state.log_line_index += 1;
        end
        close(file);
    end
end

function log_dump_prob(state::FinchState)
    if state.use_log && state.config.proc_rank == 0
        file = open(state.log_file, "a");
        println(file, string(state.log_line_index)*".\tDumping problem details:");
        state.log_line_index += 1;
        for f in fieldnames(FinchProblem)
            println(file, string(state.log_line_index)*".\t\t"*string(f)*" = "*string(getfield(p, f)));
            state.log_line_index += 1;
        end
        close(file);
    end
end

function close_log(state::FinchState)
    if state.use_log && state.config.proc_rank == 0
        log_entry(state, "Completed. Closing Log.", 0);
        state.use_log = false;
    end
end