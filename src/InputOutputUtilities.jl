# turn a codon set into a string representation for printing
function codon_set_to_str(codon_set::Vector{LongDNA{4}})
    # validate codon set
    _validate_codon_set(codon_set)

    codon_str = String.(codon_set)
    formatted_str = "\"" * join(codon_str, "\", \"") * "\""
    return formatted_str
end


# parse one compact CSV line "COD1|COD2,idx1|idx2" to Vector{LongDNA{4}}
function get_codon_set_from_line(line::AbstractString)
    # do not allow empty line
    isempty(line) && throw(ArgumentError("line is empty."))

    # require exactly one comma separator
    count(==(','), line) == 1 || throw(
        ArgumentError("line has wrong format. Expected format: \"COD1|COD2|...|CODn,idx1|idx2|...|idxn\"."),
    )
    parts = split(line, ","; limit = 2)

    # check codon list format
    codon_tokens = split(strip(parts[1]), "|"; keepempty = false)
    # turn codon_tokens into codon_set
    codon_set = LongDNA{4}.(codon_tokens)
    # validate codon_set
    _validate_codon_set(codon_set)

    # check index list format
    idx_tokens = split(strip(parts[2]), "|"; keepempty = false)
    # do not allow empty index list
    isempty(idx_tokens) && throw(ArgumentError("index list is missing."))
    # indices must be positive integers
    all(t -> occursin(r"^\d+$", t), idx_tokens) ||
        throw(ArgumentError("index list must contain only positive integers."))

    # check that codon and index counts match
    length(idx_tokens) == length(codon_tokens) || throw(
        ArgumentError("codon and index counts differ ($(length(codon_tokens)) vs $(length(idx_tokens)))."),
    )

    # check that idx_tokens are equivalent to combination indices
    comb_idxs_line = parse.(Int, idx_tokens)
    comb_idxs_codon_set = _get_comb_from_codon_set(codon_set)
    comb_idxs_line == comb_idxs_codon_set || throw(
        ArgumentError(
            "index list does not match codon list. Expected indices: $(comb_idxs_codon_set), got: $(comb_idxs_line).",
        ),
    )

    return codon_set
end


# write codon set and combination indices to CSV line: "COD1|COD2,idx1|idx2"
function write_res(io::IO, codon_set::Vector{LongDNA{4}}, comb::Vector{Int})
    # validate codon set
    _validate_codon_set(codon_set)
    # validate combination indices
    _validate_comb(comb)

    # compact CSV-style line: "COD1|COD2,idx1|idx2"
    codon_str = string.(codon_set)
    idx_str = string.(comb)
    println(io, join(codon_str, "|"), ",", join(idx_str, "|"))
    return true
end


# get corresponding combination from codon set
function _get_comb_from_codon_set(codon_set::Vector{LongDNA{4}})
    # validate codon set
    _validate_codon_set(codon_set)

    # get indices of codons in codon_set
    idxs = getindex.(Ref(CODON_INDEX), codon_set)
    return idxs
end
