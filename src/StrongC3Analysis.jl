using Base.Threads: @threads, @spawn, ReentrantLock, nthreads, Atomic

# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# process strong c3 combinations incrementally from previous results file
function process_strong_c3_combinations_by_combination_size(
    res_dir::AbstractString,
    ckp_path::AbstractString,
    stop_flag::Atomic{Bool};
    worker_cnt::Int = nthreads(),
    debug::Bool = false,
)
    # do not allow worker count < 1
    (worker_cnt < 1) && throw(ArgumentError("Invalid input for wroker_cnt (must be greater than 0)."))

    # load checkpoint if exists
    isfile(ckp_path) || throw(ArgumentError("Checkpoint file not found: $ckp_path"))
    checkpoint = _load_ckp(ckp_path)
    comb_size = checkpoint.comb_size
    # do not allow combination_size < 1
    (comb_size < 1 || comb_size > 20) &&
        throw(ArgumentError("Invalid parameter at pos1 (must be between 1-20)."))
    curr_line = checkpoint.curr_line
    w_mode = curr_line == 0 ? "w" : "a"

    # check if result directory exists
    isdir(res_dir) || throw(ArgumentError("Invalid parameter at pos2: directory not found: $res_dir"))
    prev_res_path = joinpath(res_dir, "result_$(comb_size - 1).csv")
    res_path = joinpath(res_dir, "result_$(comb_size).csv")

    # get total lines count in previous result file
    prev_line_count = comb_size == 1 ? length(ALL_CODONS) : countlines(prev_res_path)

    try
        # special case
        if comb_size == 1
            start_idx = curr_line + 1
            open(res_path, w_mode) do io
                @inbounds for i in start_idx:prev_line_count
                    if stop_flag[]
                        println("Processing interrupted at line $curr_line for combination size $comb_size.")
                        break
                    end
                    codon = ALL_CODONS[i]
                    result_to_csv!(io, [codon], [i])
                    curr_line = i
                end
            end
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
                    for (i, line) in enumerate(eachline(input))
                        # skip lines until curr_line is reached
                        i <= curr_line && continue
                        if stop_flag[]
                            println("Processing interrupted at line $i for combination size $comb_size.")
                            break
                        end
                        prev_codon_set = extract_codon_set_from_result(line)
                        isempty(prev_codon_set) && continue
                        length(prev_codon_set) == comb_size - 1 || throw(
                            ArgumentError(
                                "Invalid codon set length in previous results file at line $i. Expected length: $(comb_size - 1), got: $(length(prev_codon_set))",
                            ),
                        )
                        put!(buffer, prev_codon_set)
                        curr_line = i
                    end
                end
                close(buffer)
            end

            open(res_path, w_mode) do output
                @sync begin
                    for _ in 1:worker_cnt
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

                                    if _is_comb_strong_c3(new_comb_idxs, debug = debug)
                                        lock(w_lock) do
                                            result_to_csv!(output, ALL_CODONS[new_comb_idxs], new_comb_idxs)
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
        end
        stop_flag[] || println("Processing strong C3 combinations of size $comb_size finished successfully.")

    catch e
        throw(e)
    finally
        if curr_line == prev_line_count
            comb_size += 1
            curr_line = 0
        end
        # save checkpoint
        _save_ckp(ckp_path, comb_size, curr_line)
        println("Checkpoint saved at path $ckp_path: combination size $comb_size, line $curr_line.")
    end
    return true
end

# convert combination to mask
@inline function _comb_to_mask(combination::Vector{Int})
    mask = UInt64(0)
    @inbounds for index in combination
        mask |= _set_codon_bit(index)
    end
    return mask
end

# build rotation masks for all codons
function _get_rot_masks(codons::Vector{LongDNA{4}})
    index_dict = Dict{LongDNA{4}, Int}()
    @inbounds for (index, codon) in enumerate(codons)
        index_dict[codon] = index
    end

    masks = Vector{UInt64}(undef, length(codons))
    @inbounds for (index, codon) in enumerate(codons)
        alpha_1 = left_shift_codon(codon, 1)
        alpha_2 = left_shift_codon(codon, 2)
        index_alpha_1 = get(index_dict, alpha_1, nothing)
        index_alpha_2 = get(index_dict, alpha_2, nothing)
        index_alpha_1 === nothing && error("First rotation not found in codons: $alpha_1 (from $codon)")
        index_alpha_2 === nothing && error("Second rotation not found in codons: $alpha_2 (from $codon)")
        masks[index] = _set_codon_bit(index_alpha_1) | _set_codon_bit(index_alpha_2)
    end
    return masks
end


# check if combination is strong C3
@inline function _is_comb_strong_c3(combination::Vector{Int}; debug::Bool = false)
    codon_set = ALL_CODONS[combination]
    data = CodonGraphData(codon_set)
    construct_graph_data!(data; debug = debug)
    return is_strong_c3(data; debug = debug)
end


# load checkpoint from CSV file (single value); returns Int or nothing
function _load_ckp(path::AbstractString)
    !isfile(path) && throw(ArgumentError("Checkpoint file not found.
    File path: $path"))
    (filesize(path) == 0) && throw(ArgumentError("Checkpoint file is empty.
    File path: $path"))
    countlines(path) != 2 &&
        throw(ArgumentError("Checkpoint file must contain exactly two lines with following format:
        comb_size,value1
        curr_line,value2
        With value1 1-20 and value2 1-n
        File path: $path"))

    pairs = Dict{String, String}()
    for line in readlines(path)
        key, value = split(line, ",")
        pairs[key] = value
    end

    # parse current_combination as Vector{Int} after removing brackets and spaces
    comb_size = parse(Int, pairs["comb_size"])
    curr_line = parse(Int, pairs["curr_line"])

    return (; comb_size, curr_line)
end


# check if mask contains any codon rotation for the given combination
@inline function _mask_has_rot(combination::Vector{Int}, mask::UInt64, rotation_masks::Vector{UInt64})
    @inbounds for index in combination
        if (mask & rotation_masks[index]) != 0
            return true
        end
    end
    return false
end


# save checkpoint as one-line CSV with column name
function _save_ckp(path::AbstractString, combination_size::Int, line_number::Int)
    mktemp() do temp_file, io
        println(io, "comb_size,", string(combination_size))
        println(io, "curr_line,", string(line_number))
        close(io)
        mv(temp_file, path; force = true)
    end
    return true
end


# set bit for a 1-based codon index
@inline _set_codon_bit(idx::Integer) = UInt64(1) << (idx - 1)
