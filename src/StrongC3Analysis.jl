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
# stream all combinations for a certain combination_size, resumable via checkpoint file
function process_strong_c3_combinations_by_combination_size(
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

    # load from checkpoint or start from scratch
    if isfile(results_path) && filesize(results_path) > 0
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
            strong_c3_count = checkpoint.strong_c3_count
            not_strong_c3_count = checkpoint.not_strong_c3_count
            processed_count = checkpoint.processed_count
            write_mode = "a"
        else
            if !isfile(checkpoint_path)
                error("Checkpoint file $checkpoint_path MISSING while results file exists and is not empty.")
            elseif filesize(checkpoint_path) == 0
                error("Checkpoint file $checkpoint_path EMPTY while results file exists and is not empty.")
            end
        end
    else
        current_combination = collect(1:combination_size)
        strong_c3_count = 0
        not_strong_c3_count = 0
        processed_count = 0
        write_mode = "w"
    end

    show_debug &&
        @debug "Starting processing combinations of size $combination_size from combination $(current_combination)..."
    try
        open(results_path, write_mode) do result_out
            if write_mode == "a"
                println(result_out, "") # ensure new line before appending
            end

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
        show_debug && @debug "Finished processing combinations of size $combination_size."
    catch err
        show_debug && @debug "Error in processing combinations of size $combination_size: $err"
        rethrow(err)
    finally
        println("FINALLY ENTERED for $combination_size.")
        # save checkpoint
        _save_strong_c3_checkpoint!(
            checkpoint_path,
            current_combination,
            processed_count,
            strong_c3_count,
            not_strong_c3_count,
        )
    end
end


# stream all combinations up to max_len, resumable via checkpoint file
function process_strong_c3_combinations(
    codons::Vector{LongDNA{4}},
    max_combination_size_length::Int,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    load_checkpoint::Bool;
    show_debug::Bool = false,
)
    # load from checkpoint or start from scratch
    check_point = load_checkpoint ? _load_strong_c3_checkpoint(checkpoint_path) : nothing
    current_combination = load_checkpoint ? check_point.current_combination : Int[1]
    processed_count = load_checkpoint ? check_point.processed_count : 0
    strong_c3_count = load_checkpoint ? check_point.strong_c3_count : 0
    not_strong_c3_count = load_checkpoint ? check_point.not_strong_c3_count : 0

    # adjust variables based on load_checkpoint
    starting_point = length(current_combination)
    write_mode = load_checkpoint ? "a" : "w"

    # get length of codon set
    length_codon_set = length(codons)

    try
        # disable default SIGINT handler to allow custom handling
        Base.exit_on_sigint(false)
        open(results_path, write_mode) do result_out
            # iterate over all combination sizes
            for combination_size in starting_point:max_combination_size_length
                while true
                    # allow interrupting the process
                    yield()

                    codon_set = codons[current_combination]

                    # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
                    if _contains_codon_rotation(codon_set)
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
        end
    catch err
        rethrow(err)
    finally
        # get last processed combination from last line in results file
        last_processed_combination = _get_last_combination_indices_from_file(results_path)
        # get next combination to be processed when resumed
        current_combination = _get_next_combination(last_processed_combination, length_codon_set)
        _save_strong_c3_checkpoint!(
            checkpoint_path,
            current_combination,
            processed_count,
            strong_c3_count,
            not_strong_c3_count,
        )
        # remove empty last line from results file
        _remove_empty_last_lines(results_path)
    end
end


# function to get 


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


# write checkpoint into path
function _save_strong_c3_checkpoint!(
    path::AbstractString,
    current_combination::Vector{Int},
    processed_count::Int,
    strong_c3_count::Int,
    not_strong_c3_count::Int,
)
    open(path, "w") do out
        println(out, "\"current_combination\" = $current_combination")
        println(
            out,
            "\"processed_count\" = $(replace(Formatting.format(processed_count; commas = true), ',' => '.'))",
        )
        println(
            out,
            "\"strong_c3_count\" = $(
                replace(Formatting.format(strong_c3_count; commas = true), ',' => '.')
            )",
        )
        println(
            out,
            "\"not_strong_c3_count\" = $(
                replace(Formatting.format(not_strong_c3_count; commas = true), ',' => '.')
            )",
        )

        println(
            out,
            "\"strong_c3_percentage\" = $(@sprintf("%.2f%%", strong_c3_count / max(processed_count, 1) * 100))",
        )
    end

    return true
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
    strong_c3_count = parse(Int, pairs["strong_c3_count"])
    not_strong_c3_count = parse(Int, pairs["not_strong_c3_count"])
    processed_count = parse(Int, pairs["processed_count"])
    strong_c3_percentage = parse(Float64, replace(pairs["strong_c3_percentage"], "%" => ""))

    # return named tuple
    return (;
        current_combination,
        strong_c3_count,
        not_strong_c3_count,
        processed_count,
        strong_c3_percentage,
    )
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
    println(idxs)
    any(isnothing, idxs) && error("Codon not found in ALL_CODONS: $(codons[findfirst(isnothing, idxs)])")

    # sort indices
    sort!(idxs)
    return idxs
end


# get next combination after current_combination
function _get_next_combination(current_combination::Vector{Int}, length_codon_set::Int)
    next_combination = copy(current_combination)
    next_combination_size = length(next_combination)
    # find the rightmost element that can be incremented
    for i in next_combination_size:-1:1
        # check if this element can be incremented
        if next_combination[i] != i + length_codon_set - next_combination_size
            next_combination[i] += 1
            # reset all elements to the right of this element
            for j in (i + 1):next_combination_size
                next_combination[j] = next_combination[j - 1] + 1
            end
            return next_combination
        end
    end

    # reached last combination of this size -> return first of next size if possible
    if next_combination_size < length_codon_set
        return collect(1:(next_combination_size + 1))
    else
        return Int[] # signal no further combinations
    end
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


# check if codon_set already contains at least two of the three cyclic rotations of any codon
function _contains_codon_rotation(codon_set::Vector{LongDNA{4}})
    for codon in codon_set
        rotation_1 = left_shift_codon(codon, 1)
        rotation_2 = left_shift_codon(codon, 2)
        # check if either rotation is in codon_lookup
        if rotation_1 in codon_set || rotation_2 in codon_set
            return true
        end
    end
    return false
end


# function to get filesize in mega bytes
function _get_filesize_mb(path::AbstractString)
    isfile(path) || error("File not found: $path")

    return filesize(path) / (1024^2)
end
