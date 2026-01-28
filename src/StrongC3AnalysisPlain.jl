using Base.Threads: Atomic
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

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
    else
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
                if _contains_rotation(codon_set)
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
        # println("FINALLY ENTERED for $combination_size.")
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


# check if codon_set already contains at least two of the three cyclic rotations of any codon
function _contains_rotation(codon_set::Vector{LongDNA{4}})
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
