# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
function csv_to_result(path)
    # Parses compact CSV-style lines: "CODON1|CODON2,...,idx1|idx2"
    results = Vector{Vector{LongDNA{4}}}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line, ","; limit = 2)
        codon_part = strip(parts[1])
        codon_strings = split(codon_part, "|"; keepempty = false)
        push!(results, LongDNA{4}.(codon_strings))
    end
    return results
end


# parse one compact CSV line "COD1|COD2,idx1|idx2" -> Vector{LongDNA{4}}
function get_codon_set_from_res(line::AbstractString)
    isempty(strip(line)) && throw(ArgumentError("Line is empty."))

    parts = split(line, ",")
    codon_set = strip(parts[1])
    codon_set_str = split(codon_set, "|"; keepempty = false)
    return LongDNA{4}.(codon_set_str)
end

# turn Vector{LongDNA{4}} into compact CSV-style string "COD1", "COD2",...
function codon_set_to_str(codon_set::Vector{LongDNA{4}})
    codon_strs = String.(codon_set)
    formatted = "\"" * join(codon_strs, "\", \"") * "\""
    return formatted
end


# write codon set and combination indices to CSV line: "COD1|COD2,idx1|idx2"
function result_to_csv!(io::IOStream, codon_set::Vector{LongDNA{4}}, comb_idxs::Vector{Int})
    # compact CSV-style line: "COD1|COD2,idx1|idx2"
    codon_str = string.(codon_set)
    idx_str = string.(comb_idxs)
    println(io, join(codon_str, "|"), ",", join(idx_str, "|"))
    return true
end
