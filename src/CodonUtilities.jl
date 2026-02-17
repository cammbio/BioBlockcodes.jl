function get_comp_rev_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)
    # get complemented, reversed codon set
    complemented_reversed_codon_set = _get_comp_codon_set(_get_rev_codon_set(codon_set))

    return complemented_reversed_codon_set
end


function left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int)
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon_set
    _validate_codon_set(codon_set)

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon_set[1]))
    # shift every codon from codon_set
    shifted_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        shifted_codon = left_shift_codon(codon, shift_by)
        push!(shifted_codon_set, shifted_codon)
    end
    return shifted_codon_set
end


function left_shift_codon(codon::LongDNA{4}, shift_by::Int)
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon
    _validate_codon(codon)

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    return shifted_codon
end


#
function _get_comp_base(base::DNA)
    # do not allow bases not in ALLOWED_BASES_DNA
    !(base in ALLOWED_BASES_DNA) && throw(ArgumentError("invalid base \"$base\": must be A, C, G or T."))

    return BASE_COMPLEMENT[base]
end


#
function _get_comp_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)

    # build complemented codon set
    comp_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the complemented codon to the complemented_codons set
        push!(comp_codon_set, _get_comp_codon(codon))
    end

    return comp_codon_set
end


#
function _get_comp_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get complemented codon
    comp_codon = LongDNA{4}(_get_comp_base.(codon))
    return comp_codon
end


# 
function _get_rev_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)

    # build reversed codon set
    rev_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(rev_codon_set, _get_rev_codon(codon))
    end

    return rev_codon_set
end


#
function _get_rev_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get reversed codon
    rev_codon = codon[end:-1:1]
    return rev_codon
end
