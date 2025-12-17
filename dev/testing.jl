using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs
using BioSequences


# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
global_logger(ConsoleLogger(Logging.Info)) # deactivate


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# create test codon x0 (which is self-complementary, circular and C3)
codon_x0 =
    LongDNA{
        4,
    }.([
        "AAC",
        "AAT",
        "ACC",
        "ATC",
        "ATT",
        "CAG",
        "CTC",
        "CTG",
        "GAA",
        "GAC",
        "GAG",
        "GAT",
        "GCC",
        "GGC",
        "GGT",
        "GTA",
        "GTC",
        "GTT",
        "TAC",
        "TTC",
    ])
# first data for first graph with codon_x0
data = CodonGraphData(codon_x0; plot_title = "Codon set: $codon_x0")
construct_graph!(data; show_plot = true, show_debug = false)
# check properties of graph
is_circular(data.graph, show_debug = false)
is_comma_free(data.graph, show_debug = false)
is_self_complementary(data, show_plot = false, show_debug = false)
is_c3(data, show_plot = false, show_debug = false)



# read file per line with all 216 maximal self-complementary C3 codon_set
open("files/216_maximal_self_complementary_c3_codes.txt", "r") do f
    for (row, codon_line) in enumerate(eachline(f))
        test_data = CodonGraphData(split(codon_line); "Codon set: $(split(codon_line))")
        println("""Testing codon set #$row:
        $(test_data.codon_set)
        """)
        construct_graph!(test_data, show_plot = false, show_debug = false)
        is_circular(test_data, show_debug = false)
        is_comma_free(test_data, show_debug = false)
        is_self_complementary(test_data, show_plot = false, show_debug = false)
        is_c3(test_data, show_plot = false, show_debug = false)
        if row == 1
            break
        end
    end
end


# example from PDF
example_codon_set = LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"])
# example for cycle detection
example_data = CodonGraphData(example_codon_set)

reverse_data = CodonGraphData(get_reverse_codon_set(example_codon_set))

alpha_1_data = CodonGraphData(left_shift_codon_set(example_codon_set, 1))

alpha_2_data = CodonGraphData(left_shift_codon_set(example_codon_set, 2))

manually_adjusted_data = CodonGraphData(example_codon_set)

construct_graph!(example_data; show_plot = true, show_debug = false)
construct_graph!(reverse_data; show_plot = true, show_debug = false)
construct_graph!(alpha_1_data; show_plot = true, show_debug = false)
construct_graph!(alpha_2_data; show_plot = true, show_debug = false)
construct_graph!(manually_adjusted_data; show_plot = true, show_debug = false)

# get N₂ and N₃N₁ for each codon and add them as vertices and edges between them
for codon in data.codon_set
    n2 = codon[2:2]
    n3n1 = string(codon[3], codon[1])
    println("n2: $n2, n3n1: $n3n1")
    add_vertice_by_label!(data, n2, show_debug = true)
    add_vertice_by_label!(data, n3n1, show_debug = true)
    add_edge_by_label!(data, n2, n3n1, show_debug = true)
    add_edge_by_label!(data, n3n1, n2, show_debug = true)
end
show_graph(manually_adjusted_data; show_debug = false)

println(vcat(example_data.vertice_labels, example_data.added_vertice_labels))
println(vcat(reverse_data.vertice_labels, reverse_data.added_vertice_labels))
println(vcat(alpha_1_data.vertice_labels, alpha_1_data.added_vertice_labels))
println(vcat(alpha_2_data.vertice_labels, alpha_2_data.added_vertice_labels))
println(vcat(manually_adjusted_data.vertice_labels, manually_adjusted_data.added_vertice_labels))

data_list = [example_data, reverse_data, alpha_1_data, alpha_2_data, manually_adjusted_data]
names = ["example_data", "reverse_data", "alpha_1_data", "alpha_2_data", "manually_adjusted_data"]
show_multiple_graphs(data_list; show_debug = true)

for edge in example_data.edge_labels
    println("Edge in example_data: $(edge[1]) -> $(edge[2])")
end
for edge in reverse_data.edge_labels
    println("Edge in reverse_data: $(edge[1]) -> $(edge[2])")
end
for edge in alpha_1_data.edge_labels
    println("Edge in alpha_1_data: $(edge[1]) -> $(edge[2])")
end
for edge in alpha_2_data.edge_labels
    println("Edge in alpha_2_data: $(edge[1]) -> $(edge[2])")
end


for i in 1:(length(data_list) - 1), j in (i + 1):length(data_list)
    common_edges = intersect(data_list[i].edge_labels, data_list[j].edge_labels)
    println("""$(names[i]) compared to $(names[j])
    -> amount of common edges: $(length(common_edges))""")
    if length(common_edges) > 0
        println("list of common edges: $common_edges")
        counter = 1
        for edge in common_edges
            println("Common edge $counter: $(edge[1]) -> $(edge[2])")
            counter += 1
        end
    end
end


for name in fieldnames(typeof(manually_adjusted_data))
    println("$name => $(getfield(manually_adjusted_data, name))\n")
end


# merge two codon graphs
function merge_codon_graphs(
    data1::CodonGraphData,
    data2::CodonGraphData;
    show_debug::Bool = false,
)::CodonGraphData
    merged_data = CodonGraphData(union(data1.codon_set, data2.codon_set))
    return merged_data
end
merged_data = merge_codon_graphs(alpha_1_data, alpha_2_data; show_debug = false)
construct_graph!(merged_data; show_plot = true, show_debug = false)

show_graph(data; show_debug = false)
display_cycles(merged_data, show_debug = true)
display_cycles(data, show_debug = true)

# test data
codon_set_1 = LongDNA{4}.(["AAC", "GTT"])
test_data_1 = CodonGraphData(codon_set_1)
construct_graph!(test_data_1; show_plot = true, show_debug = false)
codon_set_2 = LongDNA{4}.(["GTA", "GTT", "GCA", "GCG"])
test_data_2 = CodonGraphData(codon_set_2)
construct_graph!(test_data_2; show_plot = true, show_debug = false)
codon_set_3 = LongDNA{4}.(["ACT", "AAA", "CGA", "CCG", "TTC", "ATA"])
test_data_3 = CodonGraphData(codon_set_3)
construct_graph!(test_data_3; show_plot = true, show_debug = false)

println(test_data_1.all_edge_labels)
println(test_data_2.all_edge_labels)
println(test_data_3.all_edge_labels)

aaa = CodonGraphData(LongDNA{4}.(["GTA", "GTT", "GCA"]))
construct_graph!(aaa; show_plot = true, show_debug = false)

for field in fieldnames(typeof(aaa))
    println("$field: $(getfield(aaa, field))\n")
end


graph = SimpleDiGraph(5)
add_edge!(graph, 1, 2)
add_edge!(graph, 2, 3)
add_edge!(graph, 4, 5)
add_edge!(graph, 5, 4)
graphplot(graph)
is_circular(graph)