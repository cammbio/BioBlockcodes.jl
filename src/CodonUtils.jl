# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
    get_complemented_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)

Return the complemented and reversed codons for every codon in `codon_set`.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Codon set to transform.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Vector{LongDNA{4}}`: Complemented and reversed codon set.

# Throws
- None.

# Example
```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
get_complemented_reversed_codon_set(codon_set)
```
"""
function get_complemented_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    temp_codon_set = get_complemented_codon_set(
        get_reversed_codon_set(codon_set; show_debug = show_debug),
        show_debug = show_debug,
    )
    show_debug && @debug "Original codon set: $(codon_set)
    -> Complemented, reversed codon set: $temp_codon_set"

    return temp_codon_set
end

"""
    get_complemented_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)

Return the complemented codons for every codon in `codon_set`.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Codon set to transform.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Vector{LongDNA{4}}`: Complemented codon set.

# Throws
- None.

# Example
```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
get_complemented_codon_set(codon_set)
```
"""
function get_complemented_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    complemented_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(complemented_codons, get_complemented_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codon_set -> complemented codon set: $complemented_codons"

    return complemented_codons
end


"""
    get_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)

Return the reversed codons for every codon in `codon_set`.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Codon set to transform.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Vector{LongDNA{4}}`: Reversed codon set.

# Throws
- None.

# Example
```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
get_reversed_codon_set(codon_set)
```
"""
function get_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
    reversed_codons = Vector{LongDNA{4}}()
    for codon in codon_set
        # add the reversed complemented codon to the reversed_codons set
        push!(reversed_codons, get_reversed_codon(codon; show_debug = show_debug))
    end
    show_debug && @debug "Original codon set: $codon_set -> reversed codon set: $reversed_codons"

    return reversed_codons
end


"""
    get_complemented_codon(codon::LongDNA{4}; show_debug::Bool = false)

Return the complemented codon. Assumes a codon length of 3.

# Arguments
- `codon::LongDNA{4}`: Codon to complement.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `LongDNA{4}`: Complemented codon.

# Throws
- None.

# Example
```julia
get_complemented_codon(LongDNA{4}("AAA"))
```
"""
function get_complemented_codon(codon::LongDNA{4}; show_debug::Bool = false)
    complemented_codon = BioSequences.complement(codon)

    if length(codon) == 3
        show_debug && @debug "Original codon: $codon, -> complemented codon: $complemented_codon"
        return complemented_codon
    end
end


"""
    get_reversed_codon(codon::LongDNA{4}; show_debug::Bool = false)

Return the reversed codon. Assumes a codon length of 3.

# Arguments
- `codon::LongDNA{4}`: Codon to reverse.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `LongDNA{4}`: Reversed codon.

# Throws
- None.

# Example
```julia
get_reversed_codon(LongDNA{4}("AAA"))
```
"""
function get_reversed_codon(codon::LongDNA{4}; show_debug::Bool = false)
    reversed_codon = reverse(codon)

    if length(codon) == 3
        show_debug && @debug "Original codon: $codon, -> reversed codon: $reversed_codon"
        return reversed_codon
    end
end


"""
    get_complemented_base(base::Char; show_debug::Bool = false)

Return the complemented base using `BASE_COMPLEMENT`.

# Arguments
- `base::Char`: Base to complement.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Char`: Complemented base.

# Throws
- `AssertionError`: If `base` is not in `BASE_COMPLEMENT`.

# Example
```julia
get_complemented_base('A')
```
"""
function get_complemented_base(base::Char; show_debug::Bool = false)
    @assert haskey(BASE_COMPLEMENT, base)
    "Base is invalid. Only A, C, G, T are allowed."
    show_debug && @debug "Original base: $base, -> complemented base: $(BASE_COMPLEMENT[base])"

    return BASE_COMPLEMENT[base]
end


"""
    left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int; show_debug::Bool = false)

Return a new codon set where each codon is left-shifted by `shift_by`.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Codon set to shift.
- `shift_by::Int`: Amount to shift.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Vector{LongDNA{4}}`: Shifted codon set.

# Throws
- None.

# Example
```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
left_shift_codon_set(codon_set, 1)
```
"""
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


"""
    left_shift_codon(codon::LongDNA{4}, shift_by::Int; show_debug::Bool = false)

Return `codon` left-shifted by `shift_by`.

# Arguments
- `codon::LongDNA{4}`: Codon to shift.
- `shift_by::Int`: Amount to shift.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `LongDNA{4}`: Shifted codon.

# Throws
- None.

# Example
```julia
left_shift_codon(LongDNA{4}("AAA"), 1)
```
"""
function left_shift_codon(codon::LongDNA{4}, shift_by::Int; show_debug::Bool = false)
    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    show_debug && @debug """Original codon: $codon
    -> shifted codon by $shift_by: $shifted_codon"""
    return shifted_codon
end
