function get_comp_rev_codon_set(codon_set::Vector{LongDNA{4}})
    # do not allow empty codon sets
    isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty."))
    # do not allow codons of length different than 3
    for codon in codon_set
        length(codon) != 3 && throw(
            ArgumentError(
                "All codons in codon set must be of length 3, got codon \"$codon\" of length $(length(codon)).",
            ),
        )
    end
    # get complemented, reversed codon set
    complemented_reversed_codon_set = _comp_codon_set(_rev_codon_set(codon_set))

    return complemented_reversed_codon_set
end


function left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int)
    # do not allow empty codon sets
    isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty."))
    # do not allow codons of length different than 3
    for codon in codon_set
        length(codon) != 3 && throw(
            ArgumentError(
                "All codons in codon set must be of length 3, got codon \"$codon\" of length $(length(codon)).",
            ),
        )
    end

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon_set[1]))
    # shift every codon from codon_set
    shifted_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        # cut of first shift_by characters and append them to the end
        shifted_codon = left_shift_codon(codon, shift_by)
        push!(shifted_codon_set, shifted_codon)
    end
    return shifted_codon_set
end


function left_shift_codon(codon::LongDNA{4}, shift_by::Int)
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    return shifted_codon
end


#
function _comp_base(base::Char)
    # do not allow bases not in BASE_COMPLEMENT
    !haskey(BASE_COMPLEMENT, base) &&
        throw(ArgumentError("Base must be one of $(keys(BASE_COMPLEMENT)), got base '$base'."))

    return BASE_COMPLEMENT[base]
end


#
function _comp_codon_set(codon_set::Vector{LongDNA{4}})
    # do not allow empty codon sets
    isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty."))
    # do not allow codons of length different than 3
    for codon in codon_set
        length(codon) != 3 && throw(
            ArgumentError(
                "All codons in codon set must be of length 3, got codon \"$codon\" of length $(length(codon)).",
            ),
        )
    end

    # build complemented codon set
    complemented_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the complemented codon to the complemented_codons set
        push!(complemented_codons, _comp_codon(codon))
    end

    return complemented_codons
end


#
function _comp_codon(codon::LongDNA{4})
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # get complemented codon
    complemented_codon = BioSequences.complement(codon)
    return complemented_codon
end


# 
function _rev_codon_set(codon_set::Vector{LongDNA{4}})
    # do not allow empty codon sets
    isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty."))
    # do not allow codons of length different than 3
    for codon in codon_set
        length(codon) != 3 && throw(
            ArgumentError(
                "All codons in codon set must be of length 3, got codon \"$codon\" of length $(length(codon)).",
            ),
        )
    end

    # build reversed codon set
    reversed_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(reversed_codons, _rev_codon(codon))
    end

    return reversed_codons
end


#
function _rev_codon(codon::LongDNA{4})
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # get reversed codon
    reversed_codon = BioSequences.reverse(codon)
    return reversed_codon
end
