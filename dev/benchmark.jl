using GCATCodes
using BioSequences: LongDNA
using Base.Threads: Atomic, nthreads, @spawn
using BenchmarkTools

# convenience aliases
const CODONS = GCATCodes.ALL_CODONS
const CANCEL = Atomic{Bool}(false)
const warmed_up = Ref(false)

function func_plain(
    codons::Vector{LongDNA{4}},
    combination_size::Int,
    cancel::Atomic{Bool};
    show_debug::Bool = false,
)
    show_debug && @debug "Entered func_plain with combination_size $combination_size"
    # get length of codon set
    length_codon_set = length(codons)
    current_combination = collect(1:combination_size)
    processed_count = 0
    strong_c3_count = 0
    not_strong_c3_count = 0


    while true
        # allow interrupting the process
        if cancel[]
            show_debug && @debug "Task with combination_size $combination_size cancelled."
            break #throw(InterruptException())
        end

        codon_set = codons[current_combination]

        # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
        if _contains_codon_rotation(codon_set)
            not_strong_c3_count += 1
        else
            # check strong C3
            data = CodonGraphData(codon_set)
            construct_graph_data!(data; show_debug = false)
            if is_strong_c3(data; show_debug = false)
                strong_c3_count += 1
            else
                not_strong_c3_count += 1
            end
        end

        processed_count += 1

        # get next combination or break if none left
        if !_increment_codon_set_combination!(current_combination, length_codon_set)
            break
        end
    end
    show_debug && @debug "Finished processing combinations of size $combination_size in func_plain."
end


function func_mask(
    codons::Vector{LongDNA{4}},
    combination_size::Int,
    cancel::Atomic{Bool};
    show_debug::Bool = false,
)
    show_debug && @debug "Entered func_mask with combination_size $combination_size"
    # get length of codon set
    length_codon_set = length(codons)
    current_combination = collect(1:combination_size)
    processed_count = 0
    strong_c3_count = 0
    not_strong_c3_count = 0

    rotation_masks = _get_rotation_masks(codons)

    while true
        # allow interrupting the process
        if cancel[]
            show_debug && @debug "Task with combination_size $combination_size cancelled."
            break #throw(InterruptException())
        end

        codon_set = codons[current_combination]
        combination_mask = _combination_to_mask(current_combination)

        # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
        if _mask_contains_codon_rotation(current_combination, combination_mask, rotation_masks)
            not_strong_c3_count += 1
        else
            # check strong C3
            data = CodonGraphData(codon_set)
            construct_graph_data!(data; show_debug = false)
            if is_strong_c3(data; show_debug = false)
                strong_c3_count += 1
            else
                not_strong_c3_count += 1
            end
        end

        processed_count += 1

        # get next combination or break if none left
        if !_increment_codon_set_combination!(current_combination, length_codon_set)
            break
        end
    end
end


function main(min_k, max_k, samples)
    println("Threads: $(nthreads()), k = $min_k:$max_k, samples = $samples")

    # warmup (JIT) auf k = 3
    if !warmed_up[]
        println("Warming up on k = 3...")
        func_plain(CODONS, 3, CANCEL; show_debug = false)
        func_mask(CODONS, 3, CANCEL; show_debug = false)
        warmed_up[] = true
    end

    for k in min_k:max_k
        println("Benchmarking k = $k...")

        trial_plain =
            @benchmark func_plain(CODONS, $k, CANCEL; show_debug = false) samples = samples evals = 1
        t_plain = median(trial_plain).time / 1e9

        trial_mask = @benchmark func_mask(CODONS, $k, CANCEL; show_debug = false) samples = samples evals = 1
        t_mask = median(trial_mask).time / 1e9


        println("plain k=$k median over $samples samples: $(t_plain) s")
        println(" mask k=$k median over $samples samples: $(t_mask) s")
        println("difference k=$k: $(t_plain - t_mask) s")
    end
    println("Benchmarking finished.")
end



function benchmark_channel_size(;
    channel_size = 1024,
    items = 2000,
    workers = Threads.nthreads(),
    work_sleep = 5e-5,
)
    ch = Channel{Int}(channel_size)
    prod_blocks = Threads.Atomic{Int}(0)

    prod_time = 0.0
    total_time = @elapsed begin
        @sync begin
            prod_task = @async begin
                t0 = time()
                for i in 1:items
                    t_put = time()
                    put!(ch, i)
                    if time() - t_put > 0
                        prod_blocks[] += 1
                    end
                end
                close(ch)
                return time() - t0
            end

            for _ in 1:workers
                @spawn begin
                    for _ in ch
                        sleep(work_sleep)
                    end
                end
            end

            prod_time = fetch(prod_task)
        end
    end

    return (
        channel_size = channel_size,
        items = items,
        workers = workers,
        work_sleep = work_sleep,
        producer_time = prod_time,
        total_time = total_time,
        prod_block_events = prod_blocks[],
    )
end


function sweep_channel_sizes(
    sizes;
    items::Int = 2000,
    workers::Int = Threads.nthreads(),
    work_sleep::Float64 = 5e-5,
)
    results = [
        benchmark_channel_size(channel_size = s, items = items, workers = workers, work_sleep = work_sleep) for s in sizes
    ]
    println("size\tprod_s\ttotal_s")
    for r in results
        println(
            "$(r.channel_size)\t$(round(r.producer_time, digits = 4))\t$(round(r.total_time, digits = 4))",
        )
    end
    return results
end

# sweep_channel_sizes([0, 64, 256, 1024, 4096]; items = 5000, workers = Threads.nthreads(), work_sleep = 5e-4)

isempty(PROGRAM_FILE) || main(1, 2, 1)

println("Done loading script.")
