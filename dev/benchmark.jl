using GCATCodes
using Base.Threads: Atomic, nthreads, @spawn
using Dates

# convenience aliases
const CODONS = GCATCodes.ALL_CODONS
const CANCEL = Atomic{Bool}(false)
const warmed_up = Ref(false)

# plain variant
function run_plain(k::Int)
    process_strong_c3_combinations_by_combination_size(
        CODONS,
        k,
        "files/bench_plain_$(k).txt",
        "files/bench_plain_$(k)_cp.txt",
        CANCEL;
        show_debug = false,
    )
end

# mask variant
function run_mask(k::Int)
    process_strong_c3_combinations_by_combination_size_with_mask(
        CODONS,
        k,
        "files/bench_mask_$(k).txt",
        "files/bench_mask_$(k)_cp.txt",
        CANCEL;
        show_debug = false,
    )
end

@inline function bench_once(fn, k::Int)
    start = time_ns()
    fn(k)
    return (time_ns() - start) / 1e9
end

function main(min_k, max_k, runs)
    println("Threads: $(nthreads()), k = $min_k:$max_k")

    # warmup (JIT) auf k = 3
    if !warmed_up[]
        println("Warming up on k = 3 ...")
        run_plain(3)
        run_mask(3)
        warmed_up[] = true
    end

    tasks = Task[]
    for k in min_k:max_k
        # delete previous results and checkpoints
        for path in (
            "files/bench_plain_$(k).txt",
            "files/bench_plain_$(k)_cp.txt",
            "files/bench_mask_$(k).txt",
            "files/bench_mask_$(k)_cp.txt",
        )
            isfile(path) && rm(path)
        end

        push!(tasks, @spawn begin
            t = bench_once(run_plain, k)
            println("plain k=$k finished in $(round(t; digits = 3)) s")
        end)

        push!(tasks, @spawn begin
            t = bench_once(run_mask, k)
            println("mask  k=$k finished in $(round(t; digits = 3)) s")
        end)
    end

    fetch.(tasks)
end

isempty(PROGRAM_FILE) || main(1, 2)
