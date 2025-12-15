# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------
const BASE_COMPLEMENT = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# returns the reversed complement codon set
function get_complement_reversed_codon_set(data::CodonGraphData; show_debug::Bool = false)
    temp_codon_set = get_complement_codon_set(
        get_reverse_codon_set(data.codon_set; show_debug = show_debug),
        show_debug = show_debug,
    )
    show_debug && @debug "Original codon set: $(data.codon_set)
    -> Complemented, reversed codon set: $temp_codon_set"

    return temp_codon_set
end

# returns the complement codon set
function get_complement_codon_set(codons::Vector{String}; show_debug::Bool = false)
    complement_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(complement_codons, get_complement_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codons -> complement codon set: $complement_codons"

    return complement_codons
end


# returns the reversed codon set
function get_reverse_codon_set(codons::Vector{String}; show_debug::Bool = false)
    reversed_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(reversed_codons, get_reverse_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codons -> reversed codon set: $reversed_codons"

    return reversed_codons
end


# returns the complement codon
function get_complement_codon(codon::AbstractString; show_debug::Bool = false)
    @assert all(haskey(BASE_COMPLEMENT, c) for c in codon)
    "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        show_debug &&
            @debug "Original codon: $codon, -> complement codon: $(String([get_complement_base(c) for c in codon]))"

        return String([get_complement_base(c) for c in codon])
    end
end


# returns the reversed codon
function get_reverse_codon(codon::AbstractString; show_debug::Bool = false)
    @assert all(haskey(BASE_COMPLEMENT, c) for c in codon)
    "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        reversed = reverse(codon)
        show_debug && @debug "Original codon: $codon, -> reversed codon: $reversed"

        return reversed
    end
end


# returns the complement base
function get_complement_base(base::Char; show_debug::Bool = false)
    @assert haskey(BASE_COMPLEMENT, base)
    "Base is invalid. Only A, C, G, T are allowed."
    show_debug && @debug "Original base: $base, -> complement base: $(BASE_COMPLEMENT[base])"

    return BASE_COMPLEMENT[base]
end