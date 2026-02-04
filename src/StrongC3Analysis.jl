using Base.Threads: @threads, @spawn, ReentrantLock, nthreads, Atomic

# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# process strong c3 combinations incrementally from previous results file
function process_strong_c3_combinations_by_combination_size(
    combination_size::Int,
    previous_results_path::AbstractString,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    stop_flag::Atomic{Bool};
    worker_count::Int = nthreads(),
    show_debug::Bool = false,
)
    # do not allow combination_size < 1
    combination_size < 1 && throw(ArgumentError("combination_size must be at least 1."))
    # do not allow combination_size > length(ALL_CODONS)
    combination_size > length(ALL_CODONS) &&
        throw(ArgumentError("combination_size must not exceed $(length(ALL_CODONS))."))

    # special case
    if combination_size == 1
        open(results_path, "w") do io
            for codon in ALL_CODONS
                # write every codon as single-entry strong C3 set
                result_to_csv!(io, [codon], [findfirst(==(codon), ALL_CODONS)])
            end
        end
    else
        # check previous results file exists
        isfile(previous_results_path) ||
            throw(ArgumentError("Previous result file not found: $previous_results_path"))

        rotation_masks = _get_rotation_masks(ALL_CODONS)
        index_lookup = Dict{LongDNA{4}, Int}(codon => idx for (idx, codon) in enumerate(ALL_CODONS))
        write_lock = ReentrantLock()

        start_line = _load_checkpoint(checkpoint_path)
        start_line === nothing && (start_line = 1)
        last_read_line = Base.Threads.Atomic{Int}(start_line - 1)

        # producer task to read previous results file into a channel
        buffer = Channel{Tuple{Int, Vector{LongDNA{4}}}}(1024)
        producer = @spawn begin
            line_counter = 0
            open(previous_results_path, "r") do input
                for line in eachline(input)
                    stop_flag[] && break
                    line_counter += 1
                    line_counter < start_line && continue
                    codon_set = extract_codon_set_from_result(line)
                    isempty(codon_set) && continue
                    length(codon_set) == combination_size - 1 || continue
                    last_read_line[] = line_counter
                    put!(buffer, (line_counter, codon_set))
                end
            end
            close(buffer)
        end

        open(results_path, "w") do output
            @sync begin
                for _ in 1:worker_count
                    @spawn begin
                        for (_line_id, base_combo) in buffer
                            base_indices = sort(getindex.(Ref(index_lookup), base_combo))
                            last_idx = base_indices[end]

                            @inbounds for new_idx in (last_idx + 1):length(ALL_CODONS)
                                combo_indices = Vector{Int}(undef, combination_size)
                                for (i, v) in enumerate(base_indices)
                                    combo_indices[i] = v
                                end
                                combo_indices[end] = new_idx
                                combo_mask = _combination_to_mask(combo_indices)
                                _mask_contains_rotation(combo_indices, combo_mask, rotation_masks) && continue

                                if _is_combination_strong_c3(combo_indices, show_debug = show_debug)
                                    lock(write_lock) do
                                        result_to_csv!(output, ALL_CODONS[combo_indices], combo_indices)
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
        if stop_flag[]
            _save_checkpoint(checkpoint_path, last_read_line[] + 1)
        else
            isfile(checkpoint_path) && rm(checkpoint_path; force = true)
        end
    end
    return true
end


# check if combination is strong C3
@inline function _is_combination_strong_c3(combination::Vector{Int}; show_debug::Bool = false)
    codon_set = ALL_CODONS[combination]
    data = CodonGraphData(codon_set)
    construct_graph_data!(data; show_debug = show_debug)
    return is_strong_c3(data; show_debug = show_debug)
end


# load checkpoint from CSV file (single value); returns Int or nothing
function _load_checkpoint(path::AbstractString)
    !isfile(path) && throw(ArgumentError("Checkpoint file not found: $path"))
    (filesize(path) == 0) && throw(ArgumentError("Checkpoint file is empty: $path"))
    countlines(path) != 1 && throw(ArgumentError("Checkpoint file must contain exactly one line: $path"))

    raw = strip(readline(path))
    value_str = occursin(',', raw) ? strip(last(split(raw, ','))) : raw
    try
        return parse(Int, value_str)
    catch
        @warn "Checkpoint file invalid, ignoring: $path"
        return nothing
    end
end


# save checkpoint as one-line CSV with column name
function _save_checkpoint(path::AbstractString, next_line::Int)
    mktemp() do tmp, io
        write(io, "line_number,", string(next_line), "\n")
        close(io)
        mv(tmp, path; force = true)
    end
    return nothing
end
