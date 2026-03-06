"""
    calc_strong_c3_comb_by_size(comb_size::Int, ckp_path::AbstractString, prev_res_path::AbstractString, res_path::AbstractString, stop_flag::Base.Threads.Atomic{Bool}; worker_count::Int=nthreads()) -> Bool

Computes strong-C3 combinations of a given size with checkpoint and resume logic.

# Arguments

  - `comb_size::Int`: Size of combinations to compute (1 to 20).
  - `ckp_path::AbstractString`: Path to the checkpoint file.
  - `prev_res_path::AbstractString`: Path to the result file of the previous combination size.
  - `res_path::AbstractString`: Path to the output file for current results.
  - `stop_flag::Base.Threads.Atomic{Bool}`: Atomic stop flag for graceful interruption.

# Keyword Arguments

  - `worker_count::Int=nthreads()`: Number of parallel workers for processing.

# Returns

  - `Bool`: `true` if processing ends regularly or is completed cleanly.

# Throws

  - `ArgumentError`: If inputs, checkpoint contents, or file states are invalid.

# Examples
"""
function calc_strong_c3_comb_by_size(
    comb_size::Int,
    ckp_path::AbstractString,
    prev_res_path::AbstractString,
    res_path::AbstractString,
    stop_flag::Base.Threads.Atomic{Bool};
    worker_count::Int = nthreads(),
)
    # only allow combination sizes from 1 to 20
    (comb_size < 1 || comb_size > 20) &&
        throw(ArgumentError("Invalid parameter at pos1 (must be between 1-20)."))

    # only allow worker_count > 0
    (worker_count < 1) && throw(ArgumentError("Invalid parameter worker_count (must be greater than 0)."))

    # load from checkpoint if results file exists and is not empty
    if isfile(res_path) && filesize(res_path) > 0
        # checkpoint file exists and is not empty
        if isfile(ckp_path) && filesize(ckp_path) > 0
            # check previous results file exists and is not empty if comb_size != 1
            if comb_size == 1
                max_lines = length(ALL_CODONS)
            else
                isfile(prev_res_path) ||
                    throw(ArgumentError("Previous result file not found: $prev_res_path"))
                filesize(prev_res_path) == 0 &&
                    throw(ArgumentError("Previous result file is empty: $prev_res_path"))
                max_lines = countlines(prev_res_path)
            end
            # load checkpoint
            ckp = _load_ckp(ckp_path)
            # validate checkpoint content
            _validate_ckp(ckp, ckp_path, comb_size, max_lines)

            next_line = ckp.next_line
            status = ckp.status

            if next_line == 0 && status == "finished"
                println("Checkpoint file indicates that combination size $comb_size is already processed.")
                return true
            end

            w_mode = "a"
        else # results file exists but checkpoint file is missing or empty - this should not happen in normal execution, but we handle it just in case
            if !isfile(ckp_path)
                throw(
                    ArgumentError(
                        "Checkpoint file \"$ckp_path\" MISSING while results file exists and not empty.",
                    ),
                )
            elseif filesize(ckp_path) == 0
                throw(
                    ArgumentError(
                        "Checkpoint file \"$ckp_path\" EMPTY while results file exists and not empty.",
                    ),
                )
            end
        end
    else # start from scratch
        next_line = 1
        status = "unfinished"
        w_mode = "w"
    end

    # flag to track if processing was cancelled during execution (e.g. due to stop flag) to provide accurate messaging and checkpoint saving
    is_cancelled = false
    # special case comb_size == 1: no need to read from previous results file, just write all codons as single combinations
    if comb_size == 1
        open(res_path, w_mode) do io
            @inbounds for line_count in next_line:length(ALL_CODONS)
                # stop flag to allow graceful interruption
                if stop_flag[]
                    next_line = line_count
                    is_cancelled = true
                    break
                end

                codon = ALL_CODONS[line_count]
                _write_res(io, [codon], [line_count])
            end
        end
    else
        codon_rot_masks = _get_rot_masks(ALL_CODONS)
        idx_dict = Dict{LongDNA{4}, Int}(codon => idx for (idx, codon) in enumerate(ALL_CODONS))
        w_lock = ReentrantLock()

        # producer task to read previous results file into a channel
        buffer = Channel{Vector{LongDNA{4}}}(1024)
        producer = @spawn begin
            open(prev_res_path, "r") do input
                for (line_count, line) in enumerate(eachline(input))
                    # skip lines until next_line is reached
                    line_count < next_line && continue

                    # stop flag to allow graceful interruption
                    if stop_flag[]
                        next_line = line_count
                        is_cancelled = true
                        break
                    end

                    prev_codon_set = get_codon_set_from_line(line)
                    # validate codon set length in previous results file
                    length(prev_codon_set) == comb_size - 1 || throw(
                        ArgumentError(
                            "invalid codon set length in previous results file at line $line_count. Expected length: $(comb_size - 1), got: $(length(prev_codon_set))",
                        ),
                    )
                    put!(buffer, prev_codon_set)
                end
            end
            close(buffer)
        end

        open(res_path, w_mode) do output
            @sync begin
                for _ in 1:worker_count
                    @spawn begin
                        for prev_codon_set in buffer
                            prev_comb_idxs = getindex.(Ref(idx_dict), prev_codon_set)
                            last_idx = prev_comb_idxs[end]
                            new_comb_idxs = Vector{Int}(undef, comb_size)
                            for (prev_idx, prev_value) in enumerate(prev_comb_idxs)
                                new_comb_idxs[prev_idx] = prev_value
                            end

                            @inbounds for new_idx in (last_idx + 1):length(ALL_CODONS)
                                new_comb_idxs[end] = new_idx
                                comb_mask = _comb_to_mask(new_comb_idxs)
                                # skip if combination contains any codon rotation
                                _mask_has_rot(new_comb_idxs, comb_mask, codon_rot_masks) && continue

                                if _is_comb_strong_c3(new_comb_idxs)
                                    lock(w_lock) do
                                        _write_res(output, ALL_CODONS[new_comb_idxs], new_comb_idxs)
                                    end
                                end
                            end
                        end
                    end
                end
                # cannot finish before producer is done reading file
                wait(producer)
            end
        end
    end

    # final check if stop flag was triggered during processing
    if is_cancelled
        println(
            "Processing for combination size $comb_size was cancelled. You can resume by running the function again with the exact same parameters.",
        )
    else
        next_line = 0
        status = "finished"
        println("Processing for combination size $comb_size completed.")
    end

    # save checkpoint
    _save_ckp(ckp_path, comb_size, next_line, status)
    println(
        "Checkpoint for comb_size $comb_size saved at path $ckp_path: comb_size: $comb_size, next_line: $next_line, status: $status.",
    )
    return true
end


# convert combination to mask
@inline function _comb_to_mask(comb::Vector{Int})
    # validate combination
    _validate_comb(comb)

    mask = UInt64(0)
    @inbounds for index in comb
        mask |= _set_codon_bit(index)
    end
    return mask
end


# build rotation masks for all codons in codon_set
function _get_rot_masks(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)

    idx_dict = Dict{LongDNA{4}, Int}()
    @inbounds for (idx, codon) in enumerate(codon_set)
        idx_dict[codon] = idx
    end

    masks = Vector{UInt64}(undef, length(codon_set))
    @inbounds for (idx, codon) in enumerate(codon_set)
        alpha_1 = left_shift_codon(codon, 1)
        alpha_2 = left_shift_codon(codon, 2)
        idx_alpha_1 = get(idx_dict, alpha_1, nothing)
        idx_alpha_2 = get(idx_dict, alpha_2, nothing)
        idx_alpha_1 === nothing && error("First rotation not found in codons: $alpha_1 (from $codon)")
        idx_alpha_2 === nothing && error("Second rotation not found in codons: $alpha_2 (from $codon)")
        masks[idx] = _set_codon_bit(idx_alpha_1) | _set_codon_bit(idx_alpha_2)
    end

    return masks
end


# check if combination is strong C3
@inline function _is_comb_strong_c3(comb::Vector{Int})
    # validate combination
    _validate_comb(comb)

    codon_set = ALL_CODONS[comb]
    cgd = CodonGraphData(codon_set)
    return is_strong_c3(cgd)
end


# load checkpoint from CSV file (single value); returns Int or nothing
function _load_ckp(path::AbstractString)
    # validate checkpoint file
    _validate_ckp_file(path)

    pairs = Dict{String, String}()
    for line in readlines(path)
        key, value = split(line, ","; limit = 2)
        pairs[strip(key)] = strip(value)
    end

    # parse current_combination as Vector{Int} after removing brackets and spaces
    comb_size = parse(Int, pairs["comb_size"])
    next_line = parse(Int, pairs["next_line"])
    status = pairs["status"]

    return (; comb_size, next_line, status)
end


# check if mask contains any codon rotation for the given combination
@inline function _mask_has_rot(comb::Vector{Int}, mask::UInt64, rot_masks::Vector{UInt64})
    # validate combination
    _validate_comb(comb)

    @inbounds for idx in comb
        if (mask & rot_masks[idx]) != 0
            return true
        end
    end
    return false
end


# save checkpoint as one-line CSV with column name
function _save_ckp(path::AbstractString, comb_size::Int, line_number::Int, status::AbstractString)
    # validate checkpoint directory
    _validate_dir(path)
    # check status value
    !(status in ("unfinished", "finished")) &&
        throw(ArgumentError("invalid status value: $status. Must be \"unfinished\" or \"finished\"."))

    mktemp() do temp_file, io
        println(io, "comb_size,", string(comb_size))
        println(io, "next_line,", string(line_number))
        println(io, "status,", status)
        close(io)
        mv(temp_file, path; force = true)
    end
    return true
end


# set bit for a 1-based codon index
@inline function _set_codon_bit(idx::Integer)
    # validate index
    (idx < 1 || idx > length(ALL_CODONS)) &&
        throw(ArgumentError("invalid codon index: $idx. Must be between 1 and $(length(ALL_CODONS))."))

    UInt64(1) << (idx - 1)
end


