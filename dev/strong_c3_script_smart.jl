#!/usr/bin/env julia

# Run incremental strong C3 search in parallel.
#
# Usage:
#   julia --project dev/strong_c3_script_smart.jl <combination_size> <prev_results_csv> <output_csv> [--debug]
#
# - combination_size: target size k (integer >= 2)
# - prev_results_csv: path to CSV-lines file from size k-1 run (format COD1|COD2,...,idx1|idx2)
# - output_csv      : path to write CSV-lines results of size k
# - --debug         : optional flag to enable verbose logging

using GCATCodes
using Base.Threads
using Base.Threads: Atomic

# shared cancel flag
const stop_flag = Atomic{Bool}(false)

function run_jobs(;
    min_combination_size::Int = 2,
    max_combination_size::Int = 20,
    worker_count::Int = nthreads(),
)
    stop_flag[] = false
    start_time = time()

    task = @spawn begin
        for k in min_combination_size:max_combination_size
            stop_flag[] && break
            prev_results = joinpath("files/tests/res2/result_$(k - 1).csv")
            out_results = joinpath("files/tests/res2/resultt_$(k).csv")

            if !isfile(prev_results)
                @warn "Skipping size $k: missing input file $prev_results"
                continue
            end

            println("Threads: $(nthreads()) | size=$k | input=$prev_results -> output=$out_results")
            GCATCodes.process_strong_c3_combinations_increment(
                k,
                prev_results,
                out_results;
                cancel = stop_flag,
                worker_count = worker_count,
            )
        end
        println("All requested sizes processed or cancelled.")
    end

    # wait for task to finish to measure total runtime
    try
        fetch(task)
    catch err
        err isa InterruptException && @info "Task cancelled."
        @warn "Task encountered an error: $(sprint(showerror, err))"
    end

    total = time() - start_time
    @info "run_jobs finished sizes $(min_combination_size):$(max_combination_size) in $(round(total, digits = 3)) seconds."
    return task
end

# cancel running job task created by run_jobs
function cancel_jobs(task::Task)
    stop_flag[] = true
    try
        fetch(task)
    catch err
        err isa InterruptException && @info "Task cancelled."
        @warn "Task encountered an error: $(sprint(showerror, err))"
    end
    stop_flag[] = false
end

println(
    "Script loaded. Call run_jobs(; min_combination_size=2, max_combination_size=20) and cancel_jobs(task) to stop.",
)
