using Base.Threads: @threads, @spawn, ReentrantLock, nthreads, Atomic
# process strong c3 combinations incrementally from previous results file
function calc_strong_c3_comb_by_size(
    comb_size::Int,
    ckp_path::AbstractString,
    prev_res_path::AbstractString,
    res_path::AbstractString,
    stop_flag::Atomic{Bool};
    worker_count::Int = nthreads(),
)
    # only allow combination sizes from 1 to 20
    (comb_size < 1 || comb_size > 20) &&
        throw(ArgumentError("Invalid parameter at pos1 (must be between 1-20)."))

    # only allow worker_count > 0
    (worker_count < 1) && throw(ArgumentError("Invalid parameter worker_count (must be greater than 0)."))

    # load from checkpoint if results file exists and is not empty
    if isfile(res_path) && filesize(res_path) > 0
        # validate checkpoint file exists and is not empty
        if isfile(ckp_path) && filesize(ckp_path) > 0
            # load checkpoint
            ckp = _load_ckp(ckp_path)
            # validate checkpoint
            _validate_checkpoint(ckp, ckp_path, comb_size)
            # check if checkpoint combination size matches the actual combination size
            (ckp.comb_size != comb_size) && throw(
                ArgumentError(
                    "Checkpoint file corrupted. Expected combination size: $comb_size, found: $(ckp.comb_size). Checkpoint file: $ckp_path",
                ),
            )
            curr_line = ckp.curr_line
            status = ckp.status

            if status == "finished"
                println("Checkpoint file indicates that combination size $comb_size is already processed.")
                return true
            end

            w_mode = "a"
        else # results file exists but checkpoint file is missing or empty - this should not happen if checkpoint is saved correctly after each run, so throw error to avoid silent data corruption
            if !isfile(ckp_path)
                throw(ArgumentError("Checkpoint file $ckp_path MISSING while results file exists."))
            elseif filesize(ckp_path) == 0
                throw(ArgumentError("Checkpoint file $ckp_path EMPTY while results file exists."))
            end
        end
    else # start from scratch
        curr_line = 1
        status = "unfinished"
        w_mode = "w"
    end

    # special case comb_size == 1: no need to read from previous results file, just write all codons as single combinations
    if comb_size == 1
        open(res_path, w_mode) do io
            @inbounds for line_count in curr_line:length(ALL_CODONS)
                if stop_flag[]
                    println("Processing interrupted at line $curr_line for combination size $comb_size.")
                    break
                end
                codon = ALL_CODONS[line_count]
                write_res(io, [codon], [line_count])
                curr_line = line_count
            end
        end
        status = "finished"
    else
        # check previous results file exists
        isfile(prev_res_path) || throw(ArgumentError("Previous result file not found: $prev_res_path"))

        codon_rot_masks = _get_rot_masks(ALL_CODONS)
        idx_dict = Dict{LongDNA{4}, Int}(codon => idx for (idx, codon) in enumerate(ALL_CODONS))
        w_lock = ReentrantLock()

        # producer task to read previous results file into a channel
        buffer = Channel{Vector{LongDNA{4}}}(1024)
        producer = @spawn begin
            open(prev_res_path, "r") do input
                for (line_count, line) in enumerate(eachline(input))
                    # skip lines until curr_line is reached
                    line_count < curr_line && continue
                    if stop_flag[]
                        curr_line = line_count
                        println("Processing interrupted at line $curr_line for combination size $comb_size.")
                        break
                    end
                    prev_codon_set = get_codon_set_from_line(line)
                    isempty(prev_codon_set) && continue
                    length(prev_codon_set) == comb_size - 1 || throw(
                        ArgumentError(
                            "Invalid codon set length in previous results file at line $line_count. Expected length: $(comb_size - 1), got: $(length(prev_codon_set))",
                        ),
                    )
                    put!(buffer, prev_codon_set)
                    curr_line = line_count
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
                                        write_res(output, ALL_CODONS[new_comb_idxs], new_comb_idxs)
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

    # mark as finished if not interrupted
    if !stop_flag[]
        status = "finished"
        println("Processing strong C3 combinations of size $comb_size finished successfully.")
    end

    # save checkpoint
    _save_ckp(ckp_path, comb_size, curr_line, status)
    println(
        "Checkpoint for comb_size $comb_size saved at path $ckp_path: comb_size: $comb_size, curr_line: $curr_line, status: $status.",
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
    data = CodonGraphData(codon_set)
    return is_strong_c3(data)
end


# load checkpoint from CSV file (single value); returns Int or nothing
function _load_ckp(path::AbstractString)
    # validate checkpoint file
    _validate_ckp_file(path)

    pairs = Dict{String, String}()
    for line in readlines(path)
        key, value = split(line, ",")
        pairs[key] = value
    end

    # parse current_combination as Vector{Int} after removing brackets and spaces
    comb_size = parse(Int, pairs["comb_size"])
    curr_line = parse(Int, pairs["curr_line"])
    status = pairs["status"]

    return (; comb_size, curr_line, status)
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

    mktemp() do temp_file, io
        println(io, "comb_size,", string(comb_size))
        println(io, "curr_line,", string(line_number))
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
