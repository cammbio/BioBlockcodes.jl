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
# global_logger(ConsoleLogger(Logging.Info)) # deactivate


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
show_graph(data; show_debug = false)
# check properties of graph
is_circular(data, show_debug = false)
is_comma_free(data, show_debug = false)
is_self_complementary(data, show_plot = false, show_debug = false)
is_c3(data, show_plot = false, show_debug = false)


# manually add vertices and edges to data_adjusted
data_adjusted = CodonGraphData(example_codon_set; plot_title = "Codon set: $example_codon_set")
construct_graph!(data_adjusted; show_plot = true, show_debug = false)

# get get N₂ and N₃N₁ for each codon and add them as vertices and edges between them
println("data_adjusted.codon_set: $(data.codon_set)")
for codon in data_adjusted.codon_set
    n2 = codon[2:2]
    n3n1 = string(codon[3], codon[1])
    println("n2: $n2, n3n1: $n3n1")
    add_vertice_by_label!(data_adjusted, n2, show_debug = false)
    add_vertice_by_label!(data_adjusted, n3n1, show_debug = false)
    connect_edge_by_label!(data_adjusted, n2, n3n1, show_debug = false)
end
show_graph(data_adjusted; show_debug = false)

# check properties of adjusted graph
is_circular(data_adjusted, show_debug = false)
is_comma_free(data_adjusted, show_debug = false)
is_self_complementary(data_adjusted, show_plot = false, show_debug = false)
is_c3(data_adjusted, show_plot = false, show_debug = false)




println(data.codon_set)
println(data_adjusted.codon_set)
println(data.vertice_labels)
println(data_adjusted.vertice_labels)
println("Edges in data")
for edge in edges(data.graph)
    print("$(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]), ")
end
println("Edges in data.adjusted")
for edge in edges(data_adjusted.graph)
    print(
        "$(data_adjusted.vertice_labels[src(edge)]) -> $(data_adjusted.vertice_labels[dst(edge)]), ",
    )
end











# test self complementarity
codon_x0_self_complementary = create_complement_reversed_codons(data::CodonGraphData)
codon_x0_self_complementary_sorted = sort(codon_x0_self_complementary)
println(codon_x0)
println(codon_x0_self_complementary)
println(codon_x0_self_complementary_sorted)
# second data for seconds graph with codon_x0_self_complementary OR codon_x0_self_complementary_sorted
data_self_complementary = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    codon_x0_self_complementary, # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $codon_x0_self_complementary", # plot_title
)
construct_graph!(data_self_complementary)
println("Graph is circular: ", is_circular(data))

is_self_complementary(data)

counter = 0
for edge in edges(data.graph)
    counter += 1
    println(
        "Edge $counter in graph 1: $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)])",
    )
    println(
        "Edge $counter in graph 2: $(data_self_complementary.vertice_labels[src(edge)]) -> $(data_self_complementary.vertice_labels[dst(edge)])",
    )
end


for edge in edges(data.graph)
    if has_edge_labeled(
        data_self_complementary,
        data.vertice_labels[src(edge)],
        data.vertice_labels[dst(edge)],
    )
        println(
            "Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) found in graph 2",
        )
    else
        println(
            "Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) NOT found in graph 2",
        )
    end
end


# read file per line with all 216 maximal self-complementary C3 codon_set
open("files/216_maximal_self_complementary_c3_codes.txt", "r") do f
    for (row, codon_line) in enumerate(eachline(f))
        test_data = CodonGraphData(
            Graphs.SimpleDiGraph(0), # graph
            split(codon_line), # codon_set
            Vector{String}(), # all_vertex_labels
            Vector{String}(), # base_vertex_labels
            Vector{String}(), # added_vertex_labels
            Vector{Tuple{String, String}}(), # all_edge_labels
            Vector{Tuple{String, String}}(), # base_edge_labels
            Vector{Tuple{String, String}}(), # added_edge_labels
            Dict{String, Int}(), # vertice_index
            "Codon set: $(split(codon_line))", # plot_title
        )
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
example_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    example_codon_set, # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $example_codon_set", # plot_title
)

reverse_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    get_reverse_codon_set(example_codon_set), # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $example_codon_set", # plot_title
)

alpha_1_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    left_shift_codon_set(example_codon_set, 1; show_debug = false), # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $example_codon_set", # plot_title
)

alpha_2_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    left_shift_codon_set(example_codon_set, 2; show_debug = false), # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $(left_shift_codon_set(example_codon_set, 2; show_debug = false))", # plot_title
)

manually_adjusted_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    example_codon_set, # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $example_codon_set", # plot_title
)

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
    merged_data = CodonGraphData(
        Graphs.SimpleDiGraph(0), # graph
        union(data1.codon_set, data2.codon_set), # codon_set
        Vector{String}(), # all_vertex_labels
        Vector{String}(), # base_vertex_labels
        Vector{String}(), # added_vertex_labels
        Vector{Tuple{String, String}}(), # all_edge_labels
        Vector{Tuple{String, String}}(), # base_edge_labels
        Vector{Tuple{String, String}}(), # added_edge_labels
        Dict{String, Int}(), # vertice_index
        "Codon set: $(union(data1.codon_set, data2.codon_set))", # plot_title
    )
    return merged_data
end
merged_data = merge_codon_graphs(alpha_1_data, alpha_2_data; show_debug = false)
construct_graph!(merged_data; show_plot = true, show_debug = false)

show_graph(data; show_debug = false)
display_cycles(merged_data, show_debug = true)
display_cycles(data, show_debug = true)



data3 = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    LongDNA{4}.(["ATA", "TAA"]), # codon_set
    Vector{String}(), # all_vertex_labels
    Vector{String}(), # base_vertex_labels
    Vector{String}(), # added_vertex_labels
    Vector{Tuple{String, String}}(), # all_edge_labels
    Vector{Tuple{String, String}}(), # base_edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
    "Codon set: $(LongDNA{4}.(["ATA", "TAA"]))", # plot_title
)
construct_graph!(data3; show_plot = true, show_debug = false)
is_circular(data3, show_debug = true)


c_set = LongDNA{4}.(["ATC", "GAT", "GAC"])
left_shift_codon_set(c_set, 1; show_debug = true)
@show c_set
println(c_set)
length(c_set)
allunique(c_set)


# test data
codon_set_1 = LongDNA{4}.(["AAC", "GTT"])
test_data_1 = create_data(codon_set_1)
construct_graph!(test_data_1; show_plot = true, show_debug = false)
codon_set_2 = LongDNA{4}.(["GTA", "GTT", "GCA", "GCG"])
test_data_2 = create_data(codon_set_2)
construct_graph!(test_data_2; show_plot = true, show_debug = false)
codon_set_3 = LongDNA{4}.(["ACT", "AAA", "CGA", "CCG", "TTC", "ATA"])
test_data_3 = create_data(codon_set_3)
construct_graph!(test_data_3; show_plot = true, show_debug = false)

println(test_data_1.all_edge_labels)
println(test_data_2.all_edge_labels)
println(test_data_3.all_edge_labels)

dataaa = CodonGraphData(LongDNA{4}.([""]))