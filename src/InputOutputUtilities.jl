# turn a codon set into a string representation for printing
function codon_set_to_str(codon_set::Vector{LongDNA{4}})
    codon_str = String.(codon_set)
    formatted_str = "\"" * join(codon_str, "\", \"") * "\""
    return formatted_str
end


# parse one compact CSV line "COD1|COD2,idx1|idx2" -> Vector{LongDNA{4}}
function get_codon_set_from_line(line::AbstractString)
    isempty(strip(line)) && throw(ArgumentError("Line is empty."))

    parts = split(line, ",")
    codon_set = strip(parts[1])
    codon_set_str = split(codon_set, "|"; keepempty = false)
    return LongDNA{4}.(codon_set_str)
end
