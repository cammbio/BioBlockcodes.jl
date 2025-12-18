# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# returns the reversed complemented codon set
function get_complemented_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    temp_codon_set = get_complemented_codon_set(
        get_reversed_codon_set(codon_set; show_debug = show_debug),
        show_debug = show_debug,
    )
    show_debug && @debug "Original codon set: $(codon_set)
    -> Complemented, reversed codon set: $temp_codon_set"

    return temp_codon_set
end

# returns the complemented codon set
function get_complemented_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    complemented_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(complemented_codons, get_complemented_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codon_set -> complemented codon set: $complemented_codons"

    return complemented_codons
end


# returns the reversed codon set
function get_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    reversed_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(reversed_codons, get_reversed_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codon_set -> reversed codon set: $reversed_codons"

    return reversed_codons
end


# returns the complemented codon
function get_complemented_codon(codon::LongDNA{4}; show_debug::Bool = false)
    complemented_codon = BioSequences.complement(codon)

    if length(codon) == 3
        show_debug && @debug "Original codon: $codon, -> complemented codon: $complemented_codon"
        return complemented_codon
    end
end


# returns the reversed codon
function get_reversed_codon(codon::LongDNA{4}; show_debug::Bool = false)
    reversed_codon = reverse(codon)

    if length(codon) == 3
        show_debug && @debug "Original codon: $codon, -> reversed codon: $reversed_codon"
        return reversed_codon
    end
end


# returns the complemented base
function get_complemented_base(base::Char; show_debug::Bool = false)
    @assert haskey(BASE_COMPLEMENT, base)
    "Base is invalid. Only A, C, G, T are allowed."
    show_debug && @debug "Original base: $base, -> complemented base: $(BASE_COMPLEMENT[base])"

    return BASE_COMPLEMENT[base]
end


# shift a codon set by k positions to the left
function left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int; show_debug::Bool = false)
    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon_set[1]))
    # shift every codon from codon_set
    shifted_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        # cut of first shift_by characters and append them to the end
        shifted_codon = left_shift_codon(codon, shift_by; show_debug = show_debug)
        push!(shifted_codon_set, shifted_codon)
    end
    show_debug && @debug """Original codon set: $codon_set
    -> shifted codon set by $shift_by: $shifted_codon_set"""
    return shifted_codon_set
end


# shift a codon by k positions to the left
function left_shift_codon(codon::LongDNA{4}, shift_by::Int; show_debug::Bool = false)
    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    show_debug && @debug """Original codon: $codon
    -> shifted codon by $shift_by: $shifted_codon"""
    return shifted_codon
end