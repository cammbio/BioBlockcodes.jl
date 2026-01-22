using GCATCodes
using Base.Threads

# cancel flag
const stop_flag = Atomic{Bool}(false)

# schedule jobs
function run_jobs()
    tasks = [
        @spawn begin
            # result_path = "files/results/strong_c3_cs$(k).txt"
            # checkpoint_path  = "files/checkpoints/strong_c3_cs$(k)_cp.txt"
            result_path = "files/results/test$(combination_size).txt"
            checkpoint_path = "files/checkpoints/test$(combination_size)_cp.txt"
            process_strong_c3_combinations_by_combination_size(
                GCATCodes.ALL_CODONS,
                combination_size,
                result_path,
                checkpoint_path,
                stop_flag;
                show_debug = false,
            )
        end for combination_size in 1:8
    ]
    println("Jobs scheduled.")
    return tasks
end

# cancel and reset all tasks
function cancel_and_reset!(tasks, stop_flag)
    stop_flag[] = true
    foreach(t -> istaskdone(t) || schedule(t, InterruptException()), tasks)
    for t in tasks
        try
            fetch(t)
        catch err
            err isa InterruptException || rethrow(err)
        end
    end

    # reset stop flag for future runs
    stop_flag[] = false
end

println("Ready to schedule jobs.")