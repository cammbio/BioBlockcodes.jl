"""
    get_comp_rev_codon_set(codon_set::Vector{LongDNA{4}}) -> Vector{LongDNA{4}}

Builds the complemented and reversed codon for each codon in the set.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Input codon set.

# Returns
- `Vector{LongDNA{4}}`: Complemented-and-reversed codon set in the same order.

# Throws
- `ArgumentError`: If `codon_set` is invalid.

# Examples
```jldoctest
julia> using GCATCodes

julia> get_comp_rev_codon_set([LongDNA{4}("ATG")])
1-element Vector{LongDNA{4}}:
 DNA "CAT"
```
"""
function get_comp_rev_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)
    # get complemented, reversed codon set
    comp_rev_codon_set = _get_comp_codon_set(_get_rev_codon_set(codon_set))

    return comp_rev_codon_set
end


"""
    left_shift_codon(codon::LongDNA{4}, shift_by::Int) -> LongDNA{4}

Performs a cyclic left shift on a codon.

# Arguments
- `codon::LongDNA{4}`: Codon to shift.
- `shift_by::Int`: Number of left shifts.

# Returns
- `LongDNA{4}`: Cyclically left-shifted codon.

# Throws
- `ArgumentError`: If `shift_by < 0` or `codon` is invalid.

# Examples
```jldoctest
julia> using GCATCodes

julia> left_shift_codon(LongDNA{4}("ATG"), 1) == LongDNA{4}("TGA")
true
```
"""
function left_shift_codon(codon::LongDNA{4}, shift_by::Int)
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon
    _validate_codon(codon)

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    shift_by == 0 && return copy(codon)
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    return shifted_codon
end


"""
    left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int) -> Vector{LongDNA{4}}

Performs a cyclic left shift on all codons in a set.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Codon set to shift.
- `shift_by::Int`: Number of left shifts.

# Returns
- `Vector{LongDNA{4}}`: Shifted codon set.

# Throws
- `ArgumentError`: If `shift_by < 0` or `codon_set` is invalid.

# Examples
```jldoctest
julia> using GCATCodes

julia> left_shift_codon_set([LongDNA{4}("ATG")], 1)
1-element Vector{LongDNA{4}}:
 DNA "TGA"
```
"""
function left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int)
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon_set
    _validate_codon_set(codon_set)

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon_set[1]))
    shift_by == 0 && return copy(codon_set)
    # shift every codon from codon_set
    shifted_codon_set = Vector{LongDNA{4}}(undef, length(codon_set))
    for (idx, codon) in enumerate(codon_set)
        shifted_codon_set[idx] = left_shift_codon(codon, shift_by)
    end
    return shifted_codon_set
end


function _get_comp_base(base::DNA)
    # do not allow bases not in ALLOWED_BASES_DNA
    !(base in ALLOWED_BASES_DNA) && throw(ArgumentError("invalid base \"$base\": must be A, C, G or T."))

    return BASE_COMPLEMENT[base]
end


function _get_comp_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get complemented codon
    comp_codon = LongDNA{4}(_get_comp_base.(codon))
    return comp_codon
end


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


function _get_rev_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get reversed codon
    rev_codon = codon[end:-1:1]
    return rev_codon
end


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


