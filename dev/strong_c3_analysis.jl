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

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# generate all combinations of a codon set by a specific size
function codon_combinations_per_size(codon_set::Vector{LongDNA{4}}, combination_size::Int)
    length_codon_set = length(codon_set)

    # do not allow combination_size <= 0
    combination_size <= 0 && throw(ArgumentError("combination_size cannot be <= 0"))
    # do not allow combination_size > length(codon_set)
    combination_size > length_codon_set &&
        throw(ArgumentError("combination_size is bigger than codon_set length"))

    combos = Vector{Vector{LongDNA{4}}}()
    # get first combination
    combination = collect(1:combination_size)
    push!(combos, codon_set[combination])

    # get next combinations
    while _get_next_codon_set_combination!(combination, combination_size, length_codon_set)
        push!(combos, codon_set[combination])
    end

    return combos
end


# generate the next combination of a codon set from the current combination
function _get_next_codon_set_combination!(
    combination::Vector{Int},
    combination_size::Int,
    length_codon_set::Int,
)
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


# get all c3 codon sets of a codon set collection
function get_all_c3_codon_sets(codon_set_collection::Vector{Vector{LongDNA{4}}})
    c3_codon_sets = Vector{Vector{LongDNA{4}}}()
    for codon_set in codon_set_collection
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)
        if is_c3(data; show_debug = false)
            push!(c3_codon_sets, codon_set)
        end
    end
    return c3_codon_sets
end







# check for each of the 216 codes all combinations of codons if they are strong C3
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    counter = 0
    open("files/strong_c3_codon_combinations.txt", "w") do out
        for line in eachline(f)
            counter += 1
            if counter == 5
                println("BREAK")
                break
            end
            codon_set = line_to_codon_set(line)

            redirect_stdout(out) do
                # check all combinations of codons
                println("For codon set ($counter): $codon_set: -------------------------------")
                for i in 1:5#length(codon_set)
                    codon_set_combinations = codon_combinations_per_size(codon_set, i)
                    for codon_set_combination in codon_set_combinations
                        # construct graph data
                        data = CodonGraphData(codon_set_combination; plot_title = "Codon set: $codon_set")
                        construct_graph_data!(data; show_debug = false)

                        # check if strong_c3
                        result = is_strong_c3(data; show_debug = false)
                        if result
                            codon_combination_string =
                                join("\"" .* string.(codon_set_combination) .* "\"", ", ")
                            println(
                                "Codon combination that is strong C3: $codon_combination_string ($(length(codon_set_combination)))",
                            )
                        end
                    end
                end
            end
        end
    end
    println("Codon combination analysis complete. Results written to files/strong_c3_codon_combinations.txt")
end


# find all "strong C3" codon sets TODO
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    for line in eachline(f)
        # construct graph data
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)


    end
end


# -------------------------------------------------- STRONG C3 ENUMERATION (STREAMING, RESUMABLE) --------------------------------------------------
# This section enumerates all combinations of ALL_CODONS up to length MAX_LENGTH, checks is_strong_c3,
# streams results to disk, and writes a checkpoint so you can abort/resume without redoing work.

const STRONG_C3_RESULTS_PATH = "files/strong_c3_codon_combinations.txt"
const STRONG_C3_CHECKPOINT_PATH = "files/strong_c3_checkpoint.txt"
const MAX_LENGTH::Int = 20

# read checkpoint file; returns nothing if missing
function _load_strong_c3_checkpoint(path::AbstractString)
    isfile(path) || return nothing

    pairs = Dict{String, String}()
    for line in readlines(path)
        isempty(strip(line)) && continue
        key, val = split(line, "=", limit = 2)
        pairs[strip(key)] = strip(val)
    end
    combination_size = parse(Int, get(pairs, "combination_size", "1"))
    combination_string = get(pairs, "current_combination", "")
    combination = isempty(combination_string) ? Int[] : parse.(Int, split(combination_string, ","))
    processed_count = parse(Int, get(pairs, "processed_count", "0"))
    strong_c3_count = parse(Int, get(pairs, "strong_c3_count", "0"))
    return (; combination_size, combination, processed_count, strong_c3_count)
end

# write checkpoint into path
function _save_strong_c3_checkpoint(
    path::AbstractString,
    combination_size::Int,
    current_combination::Vector{Int},
    processed_count::Int,
    strong_c3_count::Int,
)
    open(path, "w") do out
        current_combination_string = join(string.(current_combination), ",")
        println(out, "combination_size=$combination_size")
        println(out, "current_combination=$current_combination_string")
        println(out, "processed_count=$processed_count")
        println(out, "strong_c3_count=$strong_c3_count")
    end
    return true
end

# stream all combinations up to max_len, resumable via checkpoint file
function process_strong_c3_combinations(
    codons::Vector{LongDNA{4}};
    max_len::Int = 20,
    resume::Bool = true,
    save_interval::Int = 1000,
    results_path::AbstractString = STRONG_C3_RESULTS_PATH,
    checkpoint_path::AbstractString = STRONG_C3_CHECKPOINT_PATH,
    show_debug::Bool = false,
)
    # restore checkpoint if available
    if resume
        show_debug && println("Loading checkpoint... (from $checkpoint_path)")
        check_point = _load_strong_c3_checkpoint(checkpoint_path)
        starting_point = check_point.starting_point
        combination_size = check_point.combination_size
        current_combination = check_point.combination
        processed_count = check_point.processed_count
        strong_c3_count = check_point.strong_c3_count
    else
        show_debug && println("Starting from scratch...")
        starting_point = 1
        current_combination = Int[]
        processed_count = 0
        strong_c3_count = 0
    end

    # get length of codon set
    codon_set_length = length(codons)

    # append if resuming, otherwise overwrite
    write_mode = resume && check_point !== nothing ? "a" : "w"
    try
        open(results_path, write_mode) do out
            open("files/test_output.txt", "w") do test_out
                for combination_size in starting_point:max_len
                    # initialize first combination for this k
                    if isempty(current_combination)
                        current_combination = collect(1:combination_size)
                    end
                    println(test_out, "combination_size: $combination_size")
                    println(test_out, "current_combination: $current_combination")
                    println("combination_size: $combination_size")
                    println("current_combination: $current_combination")

                    while true
                        codon_set = codons[current_combination]
                        println(test_out, "current_combination: $current_combination")
                        println(test_out, "codon_set: $codon_set")
                        println("current_combination: $current_combination")
                        println("codon_set: $codon_set")

                        # check strong C3
                        data = CodonGraphData(codon_set)
                        construct_graph_data!(data; show_debug = false)
                        if is_strong_c3(data; show_debug = false)
                            codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                            println(out, "k=$combination_size strong C3: $codon_combination_string")
                            flush(out) # keep output on disk if aborted
                            strong_c3_count += 1
                        end

                        processed_count += 1

                        # advance to next current_combination; current_combination is mutated in-place
                        has_next = _get_next_codon_set_combination!(
                            current_combination,
                            combination_size,
                            codon_set_length,
                        )

                        has_next || break
                    end
                    current_combination = Int[] # reset for next k
                end
            end
        end

        return (; processed_count, strong_c3_count)
    catch err
        # on error, save checkpoint
        _save_strong_c3_checkpoint(
            checkpoint_path,
            combination_size,
            current_combination,
            processed_count,
            strong_c3_count,
        )
        rethrow(err)
    finally
    end
end

# start fresh
process_strong_c3_combinations(
    ALL_CODONS;
    max_len = 3,
    resume = false,
    save_interval = 1000,
    results_path = STRONG_C3_RESULTS_PATH,
    checkpoint_path = STRONG_C3_CHECKPOINT_PATH,
    show_debug = true,
)

# resume
process_strong_c3_combinations(
    ALL_CODONS;
    max_len = 3,
    resume = true,
    save_interval = 1000,
    results_path = STRONG_C3_RESULTS_PATH,
    checkpoint_path = STRONG_C3_CHECKPOINT_PATH,
    show_debug = true,
)