# turn a codon set into a string representation for printing
function codon_set_to_str(codon_set::Vector{LongDNA{4}})
    # do not allow empty codon set
    isempty(codon_set) && throw(ArgumentError("codon set is empty."))
    # codons must have length 3 and contain only allowed bases
    for codon in codon_set
        codon in ALL_CODONS || throw(ArgumentError("codon \"$codon\" in \"codon_set\" is not a valid codon!"))
    end

    codon_str = String.(codon_set)
    formatted_str = "\"" * join(codon_str, "\", \"") * "\""
    return formatted_str
end


# parse one compact CSV line "COD1|COD2,idx1|idx2" to Vector{LongDNA{4}}
function get_codon_set_from_line(line::AbstractString)
    # do not allow empty line
    isempty(line) && throw(ArgumentError("line is empty."))

    # do not allow more than one comma
    parts = split(line, ",")
    length(parts) == 2 || throw(
        ArgumentError("line has wrong format. Expected format: \"COD1|COD2|...|CODn,idx1|idx2|...|idxn\"."),
    )

    # check codon list format
    codon_tokens = split(strip(parts[1]), "|"; keepempty = false)
    # do not allow empty codon list
    isempty(codon_tokens) && throw(ArgumentError("codon list is empty."))
    # codons must have length 3 and contain only allowed bases
    for codon in codon_tokens
        length(codon) == 3 || throw(ArgumentError("codon \"$codon\" must have length 3."))
        all(base -> base in ALLOWED_BASES_STR, codon) ||
            throw(ArgumentError("codon \"$codon\" contains invalid base. Allowed: A, C, G, T."))
    end

    # check index list format
    idx_tokens = split(strip(parts[2]), "|"; keepempty = false)
    # do not allow empty index list
    isempty(idx_tokens) && throw(ArgumentError("Index list is missing."))
    # indices must be positive integers
    all(t -> occursin(r"^\d+$", t), idx_tokens) ||
        throw(ArgumentError("index list must contain only positive integers."))

    # check that codon and index counts match
    length(idx_tokens) == length(codon_tokens) || throw(
        ArgumentError("codon and index counts differ ($(length(codon_tokens)) vs $(length(idx_tokens)))."),
    )

    return LongDNA{4}.(codon_tokens)
end
