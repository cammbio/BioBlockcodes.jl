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


# turn file line (format: "YXX", "XXY", "YXZ") to codon set
function line_to_codon_set(line::String)
    codon_array = split(line, ", ")
    codon_array = replace.(codon_array, "\"" => "")
    codon_set = LongDNA{4}.(codon_array)
    return codon_set
end


# parse one compact CSV line "COD1|COD2,idx1|idx2" -> Vector{LongDNA{4}}
function extract_codon_set_from_result(line::AbstractString)
    isempty(strip(line)) && throw(ArgumentError("Line is empty."))

    parts = split(line, ",")
    codon_part = strip(parts[1])
    codon_strings = split(codon_part, "|"; keepempty = false)
    return codon_strings
end

# write codon set and combination indices to CSV line: "COD1|COD2,idx1|idx2"
function result_to_csv!(io, codon_set, combination_indices)
    # compact CSV-style line: "COD1|COD2,idx1|idx2"
    codon_strings = codon_set isa Vector{LongDNA{4}} ? string.(codon_set) : string.(LongDNA{4}.(codon_set))
    idx_strings = string.(combination_indices)
    write(io, join(codon_strings, "|"))
    write(io, ',')
    write(io, join(idx_strings, "|"))
    write(io, '\n')
end
