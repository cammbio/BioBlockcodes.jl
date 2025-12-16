using BioSequences
using Graphs


# -------------------------------------------------- STRUCTS --------------------------------------------------


# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    graph::Graphs.SimpleDiGraph # directed graph
    codon_set::Vector{LongDNA{4}} # codon set represented in the graph
    all_vertex_labels::Vector{String} # all vertice labels
    base_vertex_labels::Vector{String} # base vertex labels
    added_vertex_labels::Vector{String} # added vertex labels
    all_edge_labels::Vector{Tuple{String, String}} # all edge labels
    base_edge_labels::Vector{Tuple{String, String}} # base edge labels
    added_edge_labels::Vector{Tuple{String, String}} # added edge labels
    vertex_index::Dict{String, Int} # from vertice label to vertice index ("AA" => 3)
    plot_title::String # title of the plot
end
