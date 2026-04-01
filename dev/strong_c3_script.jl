using Revise
using BioBlockcodes
using Base.Threads

# stop_flag
const stop_flag = Base.Threads.Atomic{Bool}(false)


function start_process(min_comb_size::Int, max_comb_size::Int; worker_count::Int = nthreads())
    # do not allow combination sizes < 1 or > 20
    (min_comb_size < 1 || min_comb_size > 20) &&
        throw(ArgumentError("Invalid parameter min_comb_size (must be between 1-20)."))
    (max_comb_size < 1 || max_comb_size > 20) &&
        throw(ArgumentError("Invalid parameter max_comb_size (must be between 1-20)."))
    # do not allow min_comb_size > max_comb_size
    (min_comb_size > max_comb_size) &&
        throw(ArgumentError("min_comb_size cannot be greater than max_comb_size."))
    # only allow worker_count > 0
    (worker_count < 1) && throw(ArgumentError("Invalid parameter worker_count (must be greater than 0)."))

    stop_flag[] = false
    task = @spawn begin
        start_time = time()
        for comb_size in min_comb_size:max_comb_size
            stop_flag[] && break
            # prev_res_path = "files/results/res_$(comb_size - 1).csv"
            # res_path = "files/results/res_$(comb_size).csv"
            # ckp_path = "files/checkpoints/ckp_$(comb_size).csv"
            prev_res_path = "files/tests/results/res_$(comb_size - 1).csv"
            res_path = "files/tests/results/res_$(comb_size).csv"
            ckp_path = "files/tests/checkpoints/ckp_$(comb_size).csv"
            sort_path = "files/tests/results/sorted_res_$(comb_size).csv"

            println("Threads: $(nthreads()) | size=$comb_size | input=$prev_res_path -> output=$res_path")
            calc_strong_c3_comb_by_size(
                comb_size,
                ckp_path,
                prev_res_path,
                res_path,
                stop_flag;
                worker_count = worker_count,
            )
            _sort_by_indices(res_path, sort_path)
        end

        if stop_flag[]
            total_time = time() - start_time
            println("Processing interrupted by user after $(round(total_time, digits = 3)) seconds.")
        else
            total_time = time() - start_time
            println(
                "All requested sizes ($(min_comb_size)-$(max_comb_size)) processed in $(round(total_time, digits = 3)) seconds.",
            )
        end
    end
end

# cancel running job task created by run_jobs
function cancel_process(task::Task)
    stop_flag[] = true
    try
        fetch(task)
    catch err
        err isa InterruptException && @info "Task cancelled."
        @warn "Task encountered an error: $(sprint(showerror, err))"
    end
end


# temp
function _sort_by_indices(infile::AbstractString, outfile::AbstractString)
    lines = readlines(infile)

    sorted = sort(lines; lt = (a, b) -> begin
        idxs_a = parse.(Int, split(last(split(a, ',')), '|'))
        idxs_b = parse.(Int, split(last(split(b, ',')), '|'))
        idxs_a < idxs_b
    end)

    open(outfile, "w") do io
        write(io, join(sorted, "\n"))
    end

    return true
end


println("Script loaded.")
