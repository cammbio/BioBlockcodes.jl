# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# turn file line (format: "YXX", "XXY", "YXZ") to codon set
function line_to_codon_set(line::String)
    codon_array = split(line, ", ")
    codon_array = replace.(codon_array, "\"" => "")
    codon_set = LongDNA{4}.(codon_array)
    return codon_set
end


# redirect output of function to file
function print_to_file(dst_path::AbstractString, fnc::Function, args...)
    open(dst_path, "w") do input
        redirect_stdout(input) do
            fnc(args...)
        end
    end
end