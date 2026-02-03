using Base.Threads: @threads, @spawn, ReentrantLock, nthreads, Atomic

# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------

"""
    run_incremental_strong_c3_smart!(
        combination_size::Int,
        prev_results_path::AbstractString,
        results_path::AbstractString;
        all_codons::Vector{LongDNA{4}} = ALL_CODONS,
        show_debug::Bool = false,
    ) -> Bool

Read strong C3 results of size `combination_size - 1` (compact CSV-lines: COD1|COD2,...),
extend each base combination by exactly one further codon (only ascending, no duplicates)
and check the new combination with `is_strong_c3`. Positive hits are written thread-safely to
`results_path` in the same format. Runs in parallel on all available threads and can be cancelled
via `cancel::Atomic{Bool}`.
"""
function process_strong_c3_combinations_increment(
    combination_size::Int,
    prev_results_path::AbstractString,
    results_path::AbstractString;
    all_codons::Vector{LongDNA{4}} = ALL_CODONS,
    cancel::Atomic{Bool} = Base.Threads.Atomic{Bool}(false),
    show_debug::Bool = false,
    worker_count::Int = nthreads(),
    checkpoint_path::AbstractString = string(results_path, ".ckpt"),
)
    combination_size >= 2 || throw(ArgumentError("combination_size must be >= 2 (got $combination_size)"))
    isfile(prev_results_path) || throw(ArgumentError("Previous result file not found: $prev_results_path"))

    rotation_masks = _get_rotation_masks(all_codons)
    index_lookup = Dict{LongDNA{4}, Int}(codon => idx for (idx, codon) in enumerate(all_codons))
    write_lock = ReentrantLock()
    buffer = Channel{Tuple{Int, Vector{LongDNA{4}}}}(1024)

    # load checkpoint if present
    start_line = _load_stream_checkpoint(checkpoint_path)
    start_line === nothing && (start_line = 1)
    write_mode = start_line == 1 ? "w" : "a"

    # track last read line for checkpoint
    last_read = Base.Threads.Atomic{Int}(start_line - 1)

    producer = @spawn begin
        line_no = 0
        open(prev_results_path, "r") do io
            for line in eachline(io)
                cancel[] && break
                line_no += 1
                line_no < start_line && continue
                codon_set = csv_line_to_codon_set(line)
                isempty(codon_set) && continue
                length(codon_set) == combination_size - 1 || continue
                last_read[] = line_no
                put!(buffer, (line_no, codon_set))
            end
        end
        close(buffer)
    end

    open(results_path, write_mode) do io
        @sync begin
            for _ in 1:worker_count
                @spawn begin
                    for (line_id, base_combo) in buffer
                        base_indices = sort(getindex.(Ref(index_lookup), base_combo))
                        last_idx = base_indices[end]

                        @inbounds for new_idx in (last_idx + 1):length(all_codons)
                            combo_indices = Vector{Int}(undef, combination_size)
                            for (i, v) in enumerate(base_indices)
                                combo_indices[i] = v
                            end
                            combo_indices[end] = new_idx
                            combo_mask = _combination_to_mask(combo_indices)
                            _mask_contains_rotation(combo_indices, combo_mask, rotation_masks) && continue

                            if _combo_is_strong_c3!(combo_indices, all_codons; show_debug = show_debug)
                                lock(write_lock) do
                                    result_to_csv!(io, all_codons[combo_indices], combo_indices)
                                end
                            end
                        end
                    end
                end
            end
            # can't finish before producer is done reading file
            wait(producer)
        end
    end

    # save checkpoint on cancel, otherwise remove old one
    if cancel[]
        _save_stream_checkpoint!(checkpoint_path, last_read[] + 1)
    else
        isfile(checkpoint_path) && rm(checkpoint_path; force = true)
    end

    return true
end

# load next start line from checkpoint file; returns Int or nothing
function _load_stream_checkpoint(path::AbstractString)
    isfile(path) || return nothing
    for line in eachline(path)
        isempty(strip(line)) && continue
        if startswith(line, "next_line=")
            return parse(Int, split(line, "=")[2])
        end
    end
    return nothing
end

# save checkpoint with next line to process
function _save_stream_checkpoint!(path::AbstractString, next_line::Int)
    open(path, "w") do io
        println(io, "next_line=$(next_line)")
    end
    return true
end


@inline function _combo_is_strong_c3!(
    combo_indices::Vector{Int},
    codons::Vector{LongDNA{4}};
    show_debug::Bool = false,
)
    codon_set = codons[combo_indices]
    data = CodonGraphData(codon_set)
    construct_graph_data!(data; show_debug = show_debug)
    return is_strong_c3(data; show_debug = show_debug)
end
