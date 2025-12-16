# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# returns the reversed complemented codon set
function get_complemented_reversed_codon_set(data::CodonGraphData; show_debug::Bool = false)
    temp_codon_set = get_complemented_codon_set(
        get_reversed_codon_set(data.codon_set; show_debug = show_debug),
        show_debug = show_debug,
    )
    show_debug && @debug "Original codon set: $(data.codon_set)
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
    show_debug &&
        @debug "Original codon set: $codon_set -> complemented codon set: $complemented_codons"

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