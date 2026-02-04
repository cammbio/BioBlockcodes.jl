using Revise
using GCATCodes
using Base.Threads

# shared cancel flag
const cancel_flag = Base.Threads.Atomic{Bool}(false)


function start_process(min_combination_size::Int, max_combination_size::Int; worker_count::Int = nthreads())
    cancel_flag[] = false

    task = @spawn begin
        start_time = time()
        for combination_size in min_combination_size:max_combination_size
            cancel_flag[] && break
            previous_results_path = "files/results/result_$(combination_size - 1).csv"
            current_results_path = "files/results/result_$(combination_size).csv"
            checkpoint_path = "files/checkpoints/test.toml"

            if !isfile(previous_results_path)
                @warn "Skipping size $combination_size: missing input file $previous_results_path"
                continue
            end

            println(
                "Threads: $(nthreads()) | size=$combination_size | input=$previous_results_path -> output=$current_results_path",
            )
            process_strong_c3_combinations_by_combination_size(
                combination_size,
                previous_results_path,
                current_results_path,
                checkpoint_path;
                cancel_flag = cancel_flag,
                worker_count = worker_count,
            )
        end
        println("All requested sizes processed or cancelled.")
        total_time = time() - start_time
        @info "start_process finished sizes $(min_combination_size):$(max_combination_size) in $(round(total_time, digits = 3)) seconds."
    end

    # keep REPL responsive
    return Base.errormonitor(task)
end

# cancel running job task created by run_jobs
function cancel_process(task::Task)
    cancel_flag[] = true
    try
        fetch(task)
    catch err
        err isa InterruptException && @info "Task cancelled."
        @warn "Task encountered an error: $(sprint(showerror, err))"
    end
end

println("Script loaded.")