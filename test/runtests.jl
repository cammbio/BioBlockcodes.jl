test_file_list = [
    "codon_graph_analysis_test.jl",
    "codon_graph_plotting_test.jl",
    "codon_utilities_test.jl",
    "input_output_utilities_test.jl",
    "strong_c3_test.jl",
    "types_test.jl",
    "validation_test.jl",
]

for file in test_file_list
    if isfile(file)
        println("\n\nRunning tests for $file...")
        include(file)
    else
        @warn "Test file \"$(file)\" not found. Skipping file."
    end
end
