using GCATCodes
using Test


println("Running tests...")

isfile("codon_graph_analysis_test.jl") ? include("codon_graph_analysis_test.jl") :
@warn "Test file \"codon_graph_analysis_test.jl\" not found. Skipping file."

isfile("codon_graph_plotting_test.jl") ? include("codon_graph_plotting_test.jl") :
@warn "Test file \"codon_graph_plotting_test.jl\" not found. Skipping file."

isfile("codon_utilities_test.jl") ? include("codon_utilities_test.jl") :
@warn "Test file \"codon_utilities_test.jl\" not found. Skipping file."

isfile("input_output_utilities_test.jl") ? include("input_output_utilities_test.jl") :
@warn "Test file \"input_output_utilities_test.jl\" not found. Skipping file."

isfile("strong_c3_test.jl") ? include("strong_c3_test.jl") :
@warn "Test file \"strong_c3_test.jl\" not found. Skipping file."

isfile("types_test.jl") ? include("types_test.jl") :
@warn "Test file \"types_test.jl\" not found. Skipping file."
