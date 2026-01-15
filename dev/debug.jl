using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs
using BioSequences
using NetworkLayout

codon_set_test = LongDNA{4}.(["ATA"])
data_test = CodonGraphData(codon_set_test; plot_title = "Test")
construct_graph_data!(data_test; show_debug = false)

alpha1_test = CodonGraphData(left_shift_codon_set(codon_set_test, 1); plot_title = "Alpha 1 Test")
construct_graph_data!(alpha1_test; show_debug = false)

alpha2_test = CodonGraphData(left_shift_codon_set(codon_set_test, 2); plot_title = "Alpha 2 Test")
construct_graph_data!(alpha2_test; show_debug = false)

data_list_test = [data_test, alpha1_test, alpha2_test]
# show_multiple_codon_graphs(data_list_test; show_debug = false)

is_c3(data_test; show_debug = false)
is_strong_c3(data_test; show_debug = true)