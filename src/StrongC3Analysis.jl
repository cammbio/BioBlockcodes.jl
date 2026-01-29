using Base.Threads: Atomic
using Formatting
using Printf
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------
const ALL_CODONS =
    LongDNA{
        4,
    }.([
        "AAC",
        "AAG",
        "AAT",
        "ACA",
        "ACC",
        "ACG",
        "ACT",
        "AGA",
        "AGC",
        "AGG",
        "AGT",
        "ATA",
        "ATC",
        "ATG",
        "ATT",
        "CAA",
        "CAC",
        "CAG",
        "CAT",
        "CCA",
        "CCG",
        "CCT",
        "CGA",
        "CGC",
        "CGG",
        "CGT",
        "CTA",
        "CTC",
        "CTG",
        "CTT",
        "GAA",
        "GAC",
        "GAG",
        "GAT",
        "GCA",
        "GCC",
        "GCG",
        "GCT",
        "GGA",
        "GGC",
        "GGT",
        "GTA",
        "GTC",
        "GTG",
        "GTT",
        "TAA",
        "TAC",
        "TAG",
        "TAT",
        "TCA",
        "TCC",
        "TCG",
        "TCT",
        "TGA",
        "TGC",
        "TGG",
        "TGT",
        "TTA",
        "TTC",
        "TTG",
    ])

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# process all combinations of codons of a given size with rotation mask optimization
function process_strong_c3_combinations_by_combination_size_with_mask(
    codons::Vector{LongDNA{4}},
    combination_size::Int,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    cancel::Atomic{Bool};
    show_debug::Bool = false,
)
    show_debug && @debug "Entered function with combination_size $combination_size"
    # get length of codon set
    length_codon_set = length(codons)
    rotation_masks = _get_rotation_masks(codons)

    # load from checkpoint or start from scratch
    if isfile(results_path) #&& filesize(results_path) > 0
        if isfile(checkpoint_path) && filesize(checkpoint_path) > 0
            show_debug &&
                @debug "Loading checkpoint from $checkpoint_path for combination_size $combination_size"
            checkpoint = _load_strong_c3_checkpoint(checkpoint_path)
            current_combination = checkpoint.current_combination
            if current_combination == collect(1:(combination_size + 1))
                show_debug &&
                    @debug "Checkpoint indicates all combinations of size $combination_size processed. Exiting."
                return
            end
            processed_count = checkpoint.processed_count
            strong_c3_count = checkpoint.strong_c3_count
            not_strong_c3_count = checkpoint.not_strong_c3_count
            write_mode = "a"
        else
            if !isfile(checkpoint_path)
                error("Checkpoint file $checkpoint_path MISSING while results file exists and is not empty.")
            elseif filesize(checkpoint_path) == 0
                error("Checkpoint file $checkpoint_path EMPTY while results file exists and is not empty.")
            end
        end
    else # start from scratch
        current_combination = collect(1:combination_size)
        processed_count = 0
        strong_c3_count = 0
        not_strong_c3_count = 0
        write_mode = "w"
    end

    show_debug &&
        @debug "Starting processing combinations of size $combination_size from combination $(current_combination)..."
    try
        open(results_path, write_mode) do result_out
            while true
                # yield()
                # allow interrupting the process
                if cancel[]
                    show_debug && @debug "Task with combination_size $combination_size cancelled."
                    break
                end

                codon_set = codons[current_combination]
                combination_mask = _combination_to_mask(current_combination)

                # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
                if _mask_contains_rotation(current_combination, combination_mask, rotation_masks)
                    not_strong_c3_count += 1
                else
                    # check strong C3
                    data = CodonGraphData(codon_set)
                    construct_graph_data!(data; show_debug = false)
                    if is_strong_c3(data; show_debug = false)
                        codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                        println(
                            result_out,
                            "Strong C3: $codon_combination_string with size $(length(codon_set))",
                        )
                        strong_c3_count += 1
                    else
                        not_strong_c3_count += 1
                    end
                end

                processed_count += 1

                # get next combination or break if none left
                if !_increment_codon_set_combination!(current_combination, length_codon_set)
                    current_combination = collect(1:(combination_size + 1))
                    break
                end
            end
        end
        show_debug && @debug "Successfully finished processing all combinations of size $combination_size."
    catch err
        show_debug && @debug "Error during processing combinations of size $combination_size: $err"
        rethrow(err)
    finally
        # calculate percentages
        max_combinations = binomial(length_codon_set, combination_size)
        processed_percentage = processed_count > 0 ? (processed_count / max_combinations) * 100 : 0.0
        strong_c3_percentage =
            (strong_c3_count > 0 && processed_count > 0) ? (strong_c3_count / processed_count) * 100 : 0.0
        not_strong_c3_percentage =
            (not_strong_c3_count > 0 && processed_count > 0) ? (not_strong_c3_count / processed_count) * 100 :
            0.0
        # save checkpoint
        _save_strong_c3_checkpoint!(
            checkpoint_path,
            current_combination,
            processed_count,
            strong_c3_count,
            not_strong_c3_count,
            processed_percentage,
            strong_c3_percentage,
            not_strong_c3_percentage,
        )
    end
end


# convert combination to mask
@inline function _combination_to_mask(combination::Vector{Int})
    mask = UInt64(0)
    @inbounds for index in combination
        mask |= _set_codon_bit(index)
    end
    return mask
end


# get last processed combination from results file
function _get_last_combination_indices_from_file(path::AbstractString)
    filesize(path) == 0 && throw(ArgumentError("File empty"))

    # get last line
    last_line = readlines(path)[end]

    # extract codons from last line
    codons = [m.captures[1] for m in eachmatch(r"\"([ACGT]{3})\"", last_line)]

    # find indices of codons in ALL_CODONS
    idxs = [findfirst(==(LongDNA{4}(c)), ALL_CODONS) for c in codons]
    any(isnothing, idxs) && error("Codon not found in ALL_CODONS: $(codons[findfirst(isnothing, idxs)])")

    # sort indices
    # sort!(idxs)
    return idxs
end


# get last line of a file
function _get_last_line(path)
    open(path, "r") do io
        seekend(io)
        buf = IOBuffer()
        pos = position(io)
        while pos > 0
            pos -= 1
            seek(io, pos)
            c = read(io, Char)
            c == '\n' && pos < position(io) && break
            write(buf, c)
        end
        String(reverse(String(take!(buf))))
    end
end


# format number with dots as thousands separators
function _get_number_with_commas(num::Integer)
    return replace(Formatting.format(num; commas = true), ',' => '.')
end

# build rotation masks for all codons
function _get_rotation_masks(codons::Vector{LongDNA{4}})
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


# generate the next combination of a codon set from the current combination
function _increment_codon_set_combination!(combination::Vector{Int}, length_codon_set::Int)
    combination_size = length(combination)
    # find the rightmost element that can be incremented
    for i in combination_size:-1:1
        # check if this element can be incremented
        if combination[i] != i + length_codon_set - combination_size
            combination[i] += 1
            # reset all elements to the right of this element
            for j in (i + 1):combination_size
                combination[j] = combination[j - 1] + 1
            end
            return true
        end
    end

    return false
end

# read checkpoint file; returns nothing if missing
function _load_strong_c3_checkpoint(path::AbstractString)
    isfile(path) || return nothing

    pairs = Dict{String, String}()
    for line in readlines(path)
        isempty(strip(line)) && continue
        key_raw, value_raw = split(line, "=", limit = 2)
        # remove whitespace and quotes from key
        key = strip(strip(String(key_raw)), ['"'])
        pairs[key] = strip(value_raw)
    end

    # parse current_combination as Vector{Int} after removing brackets and spaces
    current_combination =
        parse.(Int, split(replace(pairs["current_combination"], ['[', ']', ' '] => ""), ","))
    processed_count = parse(Int, replace(pairs["processed_count"], "." => ""))
    strong_c3_count = parse(Int, replace(pairs["strong_c3_count"], "." => ""))
    not_strong_c3_count = parse(Int, replace(pairs["not_strong_c3_count"], "." => ""))

    # return named tuple
    return (; current_combination, processed_count, strong_c3_count, not_strong_c3_count)
end


# check if mask contains any codon rotation for the given combination
@inline function _mask_contains_rotation(
    combination::Vector{Int},
    mask::UInt64,
    rotation_masks::Vector{UInt64},
)
    @inbounds for index in combination
        if (mask & rotation_masks[index]) != 0
            return true
        end
    end
    return false
end


# remove all empty last lines from a file
function _remove_empty_last_lines(path::AbstractString)
    open(path, "r+") do io
        seekend(io)
        pos = position(io)
        pos == 0 && return  # leere Datei

        while pos > 0
            seek(io, pos - 1)
            b = read(io, UInt8)
            if b == 0x0a # \n
                pos -= 1
                # optionales \r vor \n entfernen (CRLF)
                if pos > 0
                    seek(io, pos - 1)
                    prev = read(io, UInt8)
                    prev == 0x0d && (pos -= 1)
                end
            elseif b == 0x0d # einzelnes \r
                pos -= 1
            else
                break
            end
        end
        truncate(io, pos)
    end

    return true
end


# write checkpoint into path
function _save_strong_c3_checkpoint!(
    path::AbstractString,
    current_combination::Vector{Int},
    processed_count::Int,
    strong_c3_count::Int,
    not_strong_c3_count::Int,
    processed_percentage::Float64,
    strong_c3_percentage::Float64,
    not_strong_c3_percentage::Float64,
)
    open(path, "w") do out
        println(out, "\"current_combination\" = $current_combination")
        println(out, "\"processed_count\" = $(_get_number_with_commas(processed_count))")
        println(out, "\"strong_c3_count\" = $(_get_number_with_commas(strong_c3_count))")
        println(
            out,
            "\"not_strong_c3_count\" = $(
       _get_number_with_commas(not_strong_c3_count))",
        )
        println(out, "\"processed_percentage\" = $processed_percentage%")
        println(out, "\"strong_c3_percentage\" = $strong_c3_percentage%")
        println(out, "\"not_strong_c3_percentage\" = $not_strong_c3_percentage%")
    end

    return true
end


# set bit for a 1-based codon index
@inline _set_codon_bit(idx::Integer) = UInt64(1) << (idx - 1)