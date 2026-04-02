module BioBlockcodesGraphMakie

using CairoMakie
using GraphMakie
using NetworkLayout
using BioCodes
using BioBlockcodes

# Pull API functions into extension module
import BioBlockcodes: plot_codon_graph, plot_multiple_codon_graphs, _get_col_count
include("CodonGraphPlotting.jl") # needs Types.jl
end