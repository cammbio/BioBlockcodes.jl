# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------
BASE_COMPLEMENT = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# generate all combinations of a codon set by a specific size MAYBE DEPRECATED
function get_codon_combinations_per_size(codon_set::Vector{LongDNA{4}}, combination_size::Int)
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
    while _increment_codon_set_combination!(combination, combination_size, length_codon_set)
        push!(combos, codon_set[combination])
    end

    return combos
end


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

  - `ArgumentError`: If `codon_set` is empty.
  - `ArgumentError`: If any codon in `codon_set` is not of length 3.

# Example

```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
get_complemented_reversed_codon_set(codon_set)
```
"""
function get_complemented_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
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
    complemented_reversed_codon_set = get_complemented_codon_set(
        get_reversed_codon_set(codon_set; show_debug = show_debug),
        show_debug = show_debug,
    )
    show_debug && @debug "Original codon set: $(codon_set)
    -> Complemented, reversed codon set: $complemented_reversed_codon_set"

    return complemented_reversed_codon_set
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

  - `ArgumentError`: If `codon_set` is empty.
  - `ArgumentError`: If any codon in `codon_set` is not of length 3.

# Example

```julia
codon_set = LongDNA{4}.(["AAA", "AAC"])
get_complemented_codon_set(codon_set)
```
"""
function get_complemented_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
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

  - `ArgumentError`: If `codon_set` is empty.
  - `ArgumentError`: If any codon in `codon_set` is not of length 3.

# Example

```julia
codon_set = LongDNA{4}.(["AGT", "AAC"])
get_reversed_codon_set(codon_set)
```
"""
function get_reversed_codon_set(codon_set::Vector{LongDNA{4}}; show_debug::Bool = false)
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

  - `ArgumentError`: If `codon` is not of length 3.

# Example

```julia
get_complemented_codon(LongDNA{4}("ATC"))
```
"""
function get_complemented_codon(codon::LongDNA{4}; show_debug::Bool = false)
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # get complemented codon
    complemented_codon = BioSequences.complement(codon)
    show_debug && @debug "Original codon: $codon, -> complemented codon: $complemented_codon"
    return complemented_codon
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

  - `ArgumentError`: If `codon` is not of length 3.

# Example

```julia
get_reversed_codon(LongDNA{4}("ACA"))
```
"""
function get_reversed_codon(codon::LongDNA{4}; show_debug::Bool = false)
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # get reversed codon
    reversed_codon = BioSequences.reverse(codon)
    show_debug && @debug "Original codon: $codon, -> reversed codon: $reversed_codon"
    return reversed_codon
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

    - `ArgumentError`: If `base` is not one of the keys in `BASE_COMPLEMENT`.

# Example

```julia
get_complemented_base('A')
```
"""
function get_complemented_base(base::Char; show_debug::Bool = false)
    # do not allow bases not in BASE_COMPLEMENT
    !haskey(BASE_COMPLEMENT, base) &&
        throw(ArgumentError("Base must be one of $(keys(BASE_COMPLEMENT)), got base '$base'."))

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

  - `ArgumentError`: If `codon_set` is empty.
  - `ArgumentError`: If any codon in `codon_set` is not of length 3.

# Example

```julia
codon_set = LongDNA{4}.(["CGA", "AAC"])
left_shift_codon_set(codon_set, 1)
```
"""
function left_shift_codon_set(codon_set::Vector{LongDNA{4}}, shift_by::Int; show_debug::Bool = false)
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

  - `ArgumentError`: If `codon` is not of length 3.

# Example

```julia    # limit shift_by to length of codon
left_shift_codon(LongDNA{4}("ATC"), 1)
```
"""
function left_shift_codon(codon::LongDNA{4}, shift_by::Int; show_debug::Bool = false)
    # do not allow codons of length different than 3
    length(codon) != 3 &&
        throw(ArgumentError("Codon must be of length 3, got codon \"$codon\" of length $(length(codon))."))

    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    show_debug && @debug """Original codon: $codon
    -> shifted codon by $shift_by: $shifted_codon"""
    return shifted_codon
end