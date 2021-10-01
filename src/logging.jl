#=
# Optional log writing
=#
export printerr, init_log, set_log_level, log_entry, log_dump_config, log_dump_prob, close_log

# verbosity levels
# 0 = minimal, not very useful
# 1 = basic messages about steps to show progress
# 2 = default, details about each step
# 3 = debug, lots of info for debugging that is hard to read
log_level = 2;

function set_log_level(lev)
    global log_level = lev;
end

function printerr(msg)
    log_entry("Error: "*msg, 0);
    println("Error: "*msg);
end

function init_log(name, dir, level=2)
    global log_file = dir*"/"*name*".txt";
    global use_log = true;
    global log_level = level;
    file = open(log_file, "w");
    println(file, "######################################");
    println(file, "# Finch Log for: "*project_name);
    println(file, "######################################");
    println(file, "(verbosity = "*string(log_level)*")");
    close(file)
end

function log_entry(text, level=2)
    global log_line_index;
    global use_log;
    if use_log && (level <= log_level)
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\t"*text);
        log_line_index += 1;
        close(file);
    end
end

function log_dump_config(c=config)
    global log_line_index;
    global use_log;
    if use_log
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\tDumping configuration:");
        log_line_index += 1;
        for f in fieldnames(Finch_config)
            println(file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(c, f)));
            log_line_index += 1;
        end
        close(file);
    end
end

function log_dump_prob(p = prob)
    global log_line_index;
    global use_log;
    if use_log
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\tDumping problem details:");
        log_line_index += 1;
        for f in fieldnames(Finch_prob)
            println(file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(p, f)));
            log_line_index += 1;
        end
        close(file);
    end
end

function close_log()
    global use_log;
    if use_log
        # log_entry("Finalizing. Dumping state.", 3);
        # if log_level >= 3
        #     log_dump_config();
        #     log_dump_prob();
        # end
        
        log_entry("Completed. Closing Log.", 0);
        use_log = false;
    end
end