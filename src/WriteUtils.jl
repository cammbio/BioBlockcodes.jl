using JSON3
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


# turn result line to codon set
function result_to_codon_set(line::String)
    target_codons = r"\"([ACGT]{3})\""
    codon_set = [LongDNA{4}(m.captures[1]) for m in eachmatch(target_codons, line)]
    return codon_set
end

# parse one compact CSV line "COD1|COD2,idx1|idx2" -> Vector{LongDNA{4}}
function csv_line_to_codon_set(line::AbstractString)
    isempty(strip(line)) && return LongDNA{4}[]
    parts = split(line, ","; limit = 2)
    codon_part = strip(parts[1])
    codon_strings = split(codon_part, "|"; keepempty = false)
    return LongDNA{4}.(codon_strings)
end


# redirect output of function to file MIGHT BE DEPRECATED
function print_to_file(dst_path::AbstractString, fnc::Function, args...)
    open(dst_path, "w") do input
        redirect_stdout(input) do
            fnc(args...)
        end
    end
end


function result_to_csv!(io, codon_set, combination_indices)
    # compact CSV-style line: "COD1|COD2,idx1|idx2"
    codon_strings = codon_set isa Vector{LongDNA{4}} ? string.(codon_set) : string.(LongDNA{4}.(codon_set))
    idx_strings = string.(combination_indices)
    write(io, join(codon_strings, "|"))
    write(io, ',')
    write(io, join(idx_strings, "|"))
    write(io, '\n')
end
