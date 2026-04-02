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
julia> using BioBlockcodes

julia> codon_set = BioBlockcodes.LongDNA{4}.(["ATG", "TGA", "TCA"]);

julia> get_comp_rev_codon_set(codon_set)
3-element Vector{BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}:
 CAT
 TCA
 TGA
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
julia> using BioBlockcodes

julia> codon = BioBlockcodes.LongDNA{4}("ACT");

julia> left_shift_codon(codon, 1)
3nt DNA Sequence:
CTA
```
"""
function left_shift_codon(codon::LongDNA{4}, shift_by::Int)
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon
    _validate_codon(codon)

    return circshift(codon; k = shift_by)
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
julia> using BioBlockcodes

julia> codon_set = BioBlockcodes.LongDNA{4}.(["TTG", "TGA"]);

julia> left_shift_codon_set(codon_set, 1)
2-element Vector{BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}:
 TGT
 GAT
```
"""
function left_shift_codon_set(codon_set::Vector{T}, shift_by::Int) where {T <: LongDNA{4}}
    # do not allow negative shift_by
    shift_by < 0 && throw(ArgumentError("shift_by must be non-negative, got shift_by = $shift_by."))
    # validate codon_set
    _validate_codon_set(codon_set)
    circshift.(codon_set; k = shift_by)
end


function _get_comp_base(base::DNA)
    # do not allow bases not in ALLOWED_BASES_DNA
    !(base in ALLOWED_BASES_DNA) && throw(ArgumentError("invalid base \"$base\": must be A, C, G or T."))

    return BioSequences.complement(base)
end


function _get_comp_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get complemented codon
    return BioSequences.complement(codon)
end


function _get_comp_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)

    # build complemented codon set
    return BioSequences.complement.(codon_set)
end


function _get_rev_codon(codon::LongDNA{4})
    # validate codon
    _validate_codon(codon)

    # get reversed codon
    return BioSequences.reverse(codon)
end


function _get_rev_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon_set
    _validate_codon_set(codon_set)

    # build reversed codon set
    return BioSequences.reverse.(codon_set)
end



