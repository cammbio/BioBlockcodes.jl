using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs
using BioSequences
using NetworkLayout

# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
global_logger(ConsoleLogger(Logging.Info)) # deactivate


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# generate all combinations of a codon set by a specific size
function codon_combinations_per_size(codon_set::Vector{LongDNA{4}}, combination_size::Int)
    length_codon_set = length(codon_set)

    # do not allow combination_size <= 0
    combination_size <= 0 && throw(ArgumentError("combination_size cannot be <= 0"))
    # do not allow combination_size > length(codon_set)
    combination_size > length_codon_set &&
        throw(ArgumentError("combination_size is bigger than codon_set length"))

    combos = Vector{Vector{LongDNA{4}}}()
    # get first combination
    combination = collect(1:combination_size)
    push!(combos, codon_set[combination])

    # get next combinations
    while _get_next_codon_set_combination!(combination, combination_size, length_codon_set)
        push!(combos, codon_set[combination])
    end

    return combos
end


# generate the next combination of a codon set from the current combination
function _get_next_codon_set_combination!(
    combination::Vector{Int},
    combination_size::Int,
    length_codon_set::Int,
)
    # find the rightmost element that can be incremented
    for i in combination_size:-1:1
        # check if this element can be incremented
        if combination[i] != i + length_codon_set - combination_size
            combination[i] += 1
            # reset all elements to the right of this element
            for j in (i + 1):combination_size
                combination[j] = combination[j - 1] + 1
            end
            return true
        end
    end

    return false
end

combinations = codon_combinations_per_size(ALL_CODONS, 3)
c3_codon_sets = get_all_c3_codon_sets(possible_combinations)
counter = 0

strong_c3_counter = 0
not_strong_c3_counter = 0

for i in 1:1
    for c3_codon_set in c3_codon_sets
        counter += 1
        data = CodonGraphData(c3_codon_set; plot_title = "c3_codon_set $counter: $c3_codon_set")
        construct_graph_data!(data; show_debug = false)
        if is_strong_c3(data; show_debug = false)
            # println("Strong C3 codon set found: $c3_codon_set")
            # show_codon_graph(data; show_debug = false)
            strong_c3_counter += 1
        else
            # println("Not a strong C3 codon set: $c3_codon_set")
            not_strong_c3_counter += 1
        end
    end
    println("Strong C3 codon sets found: $strong_c3_counter")
    println("Not strong C3 codon sets found: $not_strong_c3_counter")
end

open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    strong_c3_counter = 0
    not_strong_c3_counter = 0
    for line in eachline(f)
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)
        if is_strong_c3(data; show_debug = false)
            # println("Strong C3 codon set found: $c3_codon_set")
            # show_codon_graph(data; show_debug = false)
            strong_c3_counter += 1
        else
            # println("Not a strong C3 codon set: $c3_codon_set")
            not_strong_c3_counter += 1
        end
    end
    println("Strong C3 codon sets found: $strong_c3_counter")
    println("Not strong C3 codon sets found: $not_strong_c3_counter")
end


# get all c3 codon sets of a codon set collection
function get_all_c3_codon_sets(codon_set_collection::Vector{Vector{LongDNA{4}}})
    c3_codon_sets = Vector{Vector{LongDNA{4}}}()
    for codon_set in codon_set_collection
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)
        if is_c3(data; show_debug = false)
            push!(c3_codon_sets, codon_set)
        end
    end
    return c3_codon_sets
end

# check if alpha1 and alpha2 codon graphs contain each vertice and edge from expanded codon graph
function check_alpha_1_and_alpha_2(original_data::CodonGraphData; show_debug::Bool = false)
    expanded_data = CodonGraphData(original_data.codon_set; plot_title = "Expanded graph for inclusion check")
    construct_graph_data!(expanded_data; show_debug = false)
    _expand_graph(expanded_data; show_debug = false)

    alpha1_data = CodonGraphData(
        left_shift_codon_set(original_data.codon_set, 1);
        plot_title = "Alpha 1 graph for inclusion check",
    )
    construct_graph_data!(alpha1_data; show_debug = false)

    alpha2_data = CodonGraphData(
        left_shift_codon_set(original_data.codon_set, 2);
        plot_title = "Alpha 2 graph for inclusion check",
    )
    construct_graph_data!(alpha2_data; show_debug = false)

    data_list = [expanded_data, alpha1_data, alpha2_data]
    show_multiple_codon_graphs(data_list; show_debug = false)

    # check vertices
    for vertex in vertices(expanded_data.graph)
        vertex_label = expanded_data.all_vertex_labels[vertex]
        println("Current vertex to be checked: $vertex_label")
        if haskey(alpha1_data.vertex_index, vertex_label)
            println("Vertex $vertex_label from expanded graph found in alpha1 graph.")
        end

        if haskey(alpha2_data.vertex_index, vertex_label)
            println("Vertex $vertex_label from expanded graph found in alpha2 graph.")
        end
    end
    println("-----------------------------------------------------------------------------------------------")
    println("-----------------------------------------------------------------------------------------------")
    println("-----------------------------------------------------------------------------------------------")

    # check edges after inverting arrows from expanded graph edges
    for edge in edges(expanded_data.graph)
        in_alpha_1::Bool = false
        in_alpha_2::Bool = false

        # get edges with inverted arrows
        # src_label = expanded_data.all_vertex_labels[dst(edge)]
        # dst_label = expanded_data.all_vertex_labels[src(edge)]
        src_label = expanded_data.all_vertex_labels[src(edge)]
        dst_label = expanded_data.all_vertex_labels[dst(edge)]
        println("Current edge to be checked: $src_label -> $dst_label")

        if haskey(alpha1_data.vertex_index, src_label) && haskey(alpha1_data.vertex_index, dst_label)
            src_index_alpha1 = alpha1_data.vertex_index[src_label]
            dst_index_alpha1 = alpha1_data.vertex_index[dst_label]
            if has_edge(alpha1_data.graph, src_index_alpha1, dst_index_alpha1)
                in_alpha_1 = true
            end
        end

        if haskey(alpha2_data.vertex_index, src_label) && haskey(alpha2_data.vertex_index, dst_label)
            src_index_alpha2 = alpha2_data.vertex_index[src_label]
            dst_index_alpha2 = alpha2_data.vertex_index[dst_label]
            if has_edge(alpha2_data.graph, src_index_alpha2, dst_index_alpha2)
                in_alpha_2 = true
            end
        end

        if in_alpha_1 && in_alpha_2
            # println(
            #     "Edge $src_label -> $dst_label from expanded graph found in BOTH alpha1 and alpha2 graph.",
            # )
            println("BOTH------------------------------------------------------------------------------")
        end

        if !in_alpha_1 && !in_alpha_2
            # println("Edge $src_label -> $dst_label from expanded graph NOT found in alpha1 or alpha2 graph.")
            println("NEITHER---------------------------------------------------------------------------")
        end

        if in_alpha_1 && !in_alpha_2
            # println("Edge $src_label -> $dst_label from expanded graph ONLY found in alpha1 graph.")7
            println("ALPHA 1 ONLY")
        end

        if in_alpha_2 && !in_alpha_1
            # println("Edge $src_label -> $dst_label from expanded graph ONLY found in alpha2 graph.")
            println("ALPHA 2 ONLY")
        end
    end
end

# -------------------------------------------------- TESTING --------------------------------------------------
codon_set = LongDNA{4}.(["AAC", "GTT", "AAG", "CTT", "AAT", "ATT", "ACC"])
codon_set =
    LongDNA{
        4,
    }.([
        "AAC",
        "GTT",
        "AAG",
        "CTT",
        "AAT",
        "ATT",
        "ACC",
        "GGT",
        "ACG",
        "CGT",
        "ACT",
        "AGT",
        "AGC",
        "GCT",
        "AGG",
        "CCT",
        "CCG",
        "CGG",
        "TCA",
        "TGA",
    ])
data_original = CodonGraphData(codon_set; plot_title = "Original")
construct_graph_data!(data_original; show_debug = false)
show_codon_graph(data_original; show_debug = false)
is_strong_c3(data_original; show_debug = true)
check_alpha_1_and_alpha_2(data_original; show_debug = true)

codon_set_test = LongDNA{4}.(["ATA", "GTA", "CCA", "CGC", "CGT"])
data_test = CodonGraphData(codon_set_test; plot_title = "Test")
construct_graph_data!(data_test; show_debug = false)

alpha1_test = CodonGraphData(left_shift_codon_set(codon_set_test, 1); plot_title = "Alpha 1 Test")
construct_graph_data!(alpha1_test; show_debug = false)

alpha2_test = CodonGraphData(left_shift_codon_set(codon_set_test, 2); plot_title = "Alpha 2 Test")
construct_graph_data!(alpha2_test; show_debug = false)

data_list_test = [data_test, alpha1_test, alpha2_test]
show_multiple_codon_graphs(data_list_test; show_debug = false)
show_codon_graph(data_test; show_debug = false)

is_c3(data_test; show_debug = false)
is_strong_c3(data_test; show_debug = true)

_has_cycle_longer_than(data_test.graph, 2; show_debug = false)

add_edge_by_label!(data_test, "CA", "A", show_debug = false)
add_edge_by_label!(data_test, "A", "CG", show_debug = false)
show_codon_graph(data_test; show_debug = false)

_has_cycle_longer_than(data_test.graph, 2; show_debug = false)

_expand_graph(data_test; show_debug = false)
show_codon_graph(data_test; show_debug = false)

println(get_cycles_all(data_test; show_debug = false))

copy_codon_set = copy(codon_set)
for codon in codon_set
    n2 = string(codon[2])
    n3n1 = string(codon[3], codon[1])
    # add_vertex_by_label!(data_original, n2, show_debug = false)
    # add_vertex_by_label!(data_original, n3n1, show_debug = false)
    # add_edge_by_label!(data_original, n2, n3n1, show_debug = false)
    # add_edge_by_label!(data_original, n3n1, n2, show_debug = false)

    new_codon = left_shift_codon(codon, 1)

    # check if new_codon already in codon_set
    if !(new_codon in codon_set)
        println("Adding codon $new_codon derived from $codon by left shift.")
        push!(copy_codon_set, new_codon)
    else
        println("Codon $new_codon derived from $codon by left shift already in codon set.")
    end

end
println(copy_codon_set)
println(length(copy_codon_set))

data_expanded = CodonGraphData(codon_set; plot_title = "Expanded")
construct_graph_data!(data_expanded; show_debug = false)
add_n2_n3n1_vertices_and_edges!(data_expanded; show_debug = false)

data_alpha1 = CodonGraphData(left_shift_codon_set(codon_set, 1); plot_title = "Alpha 1")
construct_graph_data!(data_alpha1; show_debug = false)

data_alpha2 = CodonGraphData(left_shift_codon_set(codon_set, 2); plot_title = "Alpha 2")
construct_graph_data!(data_alpha2; show_debug = false)

a = get_cycles_difference(data_expanded, data_original; show_debug = false)
println(a)
println(codon_combinations_per_size(codon_set, 12))
println(length(codon_combinations_per_size(codon_set, 12)))


function testing(codon_set::Vector{LongDNA{4}})
    data_original = CodonGraphData(codon_set; plot_title = "Original")
    construct_graph_data!(data_original; show_debug = false)
    # display_all_cycles(data_original; show_debug = false)
    # get_cycle_count(data_original.graph; show_debug = false)

    data_expanded = CodonGraphData(codon_set; plot_title = "Expanded")
    construct_graph_data!(data_expanded; show_debug = false)
    add_n2_n3n1_vertices_and_edges!(data_expanded; show_debug = false)
    # display_all_cycles(data_expanded; show_debug = false)
    # get_cycle_count(data_expanded.graph; show_debug = false)

    data_alpha_1 = CodonGraphData(left_shift_codon_set(codon_set, 1); plot_title = "Alpha 1")
    construct_graph_data!(data_alpha_1; show_debug = false)
    # display_all_cycles(alpha1; show_debug = false)
    # get_cycle_count(alpha1.graph; show_debug = false)

    data_alpha_2 = CodonGraphData(left_shift_codon_set(codon_set, 2); plot_title = "Alpha 2")
    construct_graph_data!(data_alpha_2; show_debug = false)
    # display_all_cycles(alpha2; show_debug = false)
    # get_cycle_count(alpha2.graph; show_debug = false)

    data_list = [data_original, data_expanded, data_alpha_1, data_alpha_2]
    show_multiple_codon_graphs(data_list; show_debug = false)
    # show_graph(data_original; show_debug = false)
    show_codon_graph(data_expanded; show_debug = false)
    # show_graph_only_cycles(data_expanded; show_debug = false)
    # show_graph(alpha1; show_debug = false)
    # show_graph(alpha2; show_debug = false)

    print_to_file("files/test_output_original.txt", get_cycles_all, data_original)
    print_to_file("files/test_output_expanded.txt", get_cycles_all, data_expanded)
    print_to_file("files/test_output_alpha1.txt", get_cycles_all, data_alpha_1)
    print_to_file("files/test_output_alpha2.txt", get_cycles_all, data_alpha_2)

end

open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    codon_set_line = readline(f)
    # for i in 1:5
    codon_set_test = line_to_codon_set(codon_set_line)[1:10]
    println("Testing codon set: $codon_set_test")
    testing(codon_set_test)
    # end
end

function show_graph_only_cycles(data::CodonGraphData; show_debug::Bool = false)
    cycles = simplecycles(data.graph)
    cycle_vertices = sort!(unique(vcat(cycles...)))
    isempty(cycle_vertices) && error("Keine Zyklen gefunden.")

    # alter -> neuer Index
    vmap = Dict(v => i for (i, v) in enumerate(cycle_vertices))
    cycle_graph = SimpleDiGraph(length(cycle_vertices))
    for cycle in cycles
        for i in eachindex(cycle)
            src = cycle[i]
            dst = cycle[mod1(i + 1, length(cycle))]
            add_edge!(cycle_graph, vmap[src], vmap[dst])
        end
    end

    graph_to_plot = cycle_graph
    labels_to_plot = [data.all_vertex_labels[v] for v in cycle_vertices]


    fig = Figure(size = (1800, 900))
    ax = Axis(
        fig[1, 1];
        xgridvisible = false,
        ygridvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
    )
    hidespines!(ax)
    ax.title = data.plot_title
    graphplot!(
        ax,
        graph_to_plot;
        layout = Spring(C = 50.0),
        nlabels = labels_to_plot,
        nlabels_color = :white,
        nlabels_size = 18,
        nlabels_offset = Point2f(0, 0),
        nlabels_align = (:center, :center),
        node_color = :black,
        node_size = 30,
        arrow_shift = :end,
        arrow_size = 12,
        edge_width = 2,
        edge_curvature = 0.9,
    )
    display(fig)
end














# -------------------------------------------------- ANALYSIS --------------------------------------------------
# group codon sets by maximal cycle length
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    maximal_cycle_length_dict = Dict{Int, Vector{Vector{LongDNA{4}}}}()
    for line in eachline(f)
        # construct graph data
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)

        # manually add N₂ and N₃N₁ vertices and connect edges
        add_n2_n3n1_vertices_and_edges!(data; show_debug = false)

        # get maximal cycle length
        maximal_cycle_length = get_max_cycle_length(data.graph; show_debug = false)
        println("Maximal cycle length: $maximal_cycle_length")
        push!(get!(maximal_cycle_length_dict, maximal_cycle_length, Vector{Vector{LongDNA{4}}}()), codon_set)
    end

    # write result to file
    open("files/codon_sets_grouped_by_maximal_cycle_length.txt", "w") do out
        redirect_stdout(out) do
            for key in sort(collect(keys(maximal_cycle_length_dict)))
                println("Maximal cycle length: $key")
                for codon_set in maximal_cycle_length_dict[key]
                    codon_set_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                    println("Codon set: $codon_set_string")
                end
                println()
            end
        end
    end
    println(
        "Codon sets grouped by maximal cycle length written to files/codon_sets_grouped_by_maximal_cycle_length.txt",
    )
end


# group codon sets by cycle count
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    cycle_count_dict = Dict{Int, Vector{Vector{LongDNA{4}}}}()
    for line in eachline(f)
        # construct graph data
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)

        # manually add N₂ and N₃N₁ vertices and connect edges
        add_n2_n3n1_vertices_and_edges!(data; show_debug = false)

        # get cycle count
        cycle_count = get_cycle_count(data.graph; show_debug = false)
        # println("Cycle count: $cycle_count")
        push!(get!(cycle_count_dict, cycle_count, Vector{Vector{LongDNA{4}}}()), codon_set)
    end

    # write result to file
    open("files/codon_sets_grouped_by_cycle_count.txt", "w") do out
        redirect_stdout(out) do
            for key in sort(collect(keys(cycle_count_dict)))
                println("Cycle count: $key, codon set count: $(length(cycle_count_dict[key]))")
                for codon_set in cycle_count_dict[key]
                    codon_set_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                    println("Codon set: $codon_set_string")
                end
                println()
            end
        end
    end
    println("Codon sets grouped by cycle count written to files/codon_sets_grouped_by_cycle_count.txt")
end


# group codon sets by cycle count for cycles of length x
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    cycle_count_length_dict = Dict{Int, Dict{Int, Vector{Vector{LongDNA{4}}}}}()
    for line in eachline(f)
        # construct graph data
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)

        # manually add N₂ and N₃N₁ vertices and connect edges
        add_n2_n3n1_vertices_and_edges!(data; show_debug = false)

        for cycle_length in 2:50
            # get cycle count of given length
            cycle_count = get_cycle_count_of_length(data, cycle_length; show_debug = false)
            if cycle_count != 0
                inner_dict_by_length =
                    get!(cycle_count_length_dict, cycle_length, Dict{Int, Vector{Vector{LongDNA{4}}}}())
                push!(get!(inner_dict_by_length, cycle_count, Vector{Vector{LongDNA{4}}}()), codon_set)
            end
        end
    end

    # write result to file
    open("files/codon_sets_grouped_by_cycle_count_and_cycle_length.txt", "w") do out
        redirect_stdout(out) do
            for length_key in sort(collect(keys(cycle_count_length_dict)))
                for count_key in sort(collect(keys(cycle_count_length_dict[length_key])))
                    println(
                        "The following $(length(cycle_count_length_dict[length_key][count_key])) codon sets contain $count_key cycles of length $length_key:",
                    )
                    for codon_set in cycle_count_length_dict[length_key][count_key]
                        codon_set_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                        println("  Codon set: $codon_set_string")
                    end
                    println()
                end
            end
        end
    end
    println(
        "Codon sets grouped by cycle count and length written to files/codon_sets_grouped_by_cycle_count_and_length.txt",
    )
end


# check for each codon in each codon set if expanding vertices and edges (per codon) keeps the graph C3
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    false_count = 0
    true_count = 0
    open("files/is_c3_after_adding_vertices_and_edges_of_only_one_codon.txt", "w") do out
        redirect_stdout(out) do
            for line in eachline(f)
                # construct graph data
                codon_set = line_to_codon_set(line)
                println("Codon set: $line")
                # manually add N₂ and N₃N₁ vertices and connect edges for each codon
                for codon in codon_set
                    data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
                    construct_graph_data!(data; show_debug = false)

                    _add_n2_n3n1_by_codon(data, codon; show_debug = false)
                    result = is_c3(data; show_debug = false)
                    if result
                        true_count += 1
                    else
                        false_count += 1
                    end
                    println("""added N₂ and N₃N₁ vertices and edges for codon: $codon
                    is_c3: $result""")
                end
            end
        end
    end
    println(
        """Codon set analysis per codon complete. Results written to files/is_c3_after_adding_vertices_and_edges_of_only_one_codon.txt
        false_count: $false_count
        true_count: $true_count""",
    )
end


# check for each of the 216 codes all combinations of codons if they are strong C3
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    counter = 0
    open("files/strong_c3_codon_combinations.txt", "w") do out
        for line in eachline(f)
            counter += 1
            if counter == 5
                println("BREAK")
                break
            end
            codon_set = line_to_codon_set(line)

            redirect_stdout(out) do
                # check all combinations of codons
                println("For codon set ($counter): $codon_set: -------------------------------")
                for i in 1:5#length(codon_set)
                    codon_set_combinations = codon_combinations_per_size(codon_set, i)
                    for codon_set_combination in codon_set_combinations
                        # construct graph data
                        data = CodonGraphData(codon_set_combination; plot_title = "Codon set: $codon_set")
                        construct_graph_data!(data; show_debug = false)

                        # check if strong_c3
                        result = is_strong_c3(data; show_debug = false)
                        if result
                            codon_combination_string =
                                join("\"" .* string.(codon_set_combination) .* "\"", ", ")
                            println(
                                "Codon combination that is strong C3: $codon_combination_string ($(length(codon_set_combination)))",
                            )
                        end
                    end
                end
            end
        end
    end
    println("Codon combination analysis complete. Results written to files/strong_c3_codon_combinations.txt")
end


# find all "strong C3" codon sets TODO
open("files/216_maximal_self_complementary_c3_codes_array.txt", "r") do f
    for line in eachline(f)
        # construct graph data
        codon_set = line_to_codon_set(line)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
        construct_graph_data!(data; show_debug = false)


    end
end



# -------------------------------------------------- INTERACTIVITY --------------------------------------------------
# graph interactivity
g = SimpleDiGraph(6)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)
add_edge!(g, 5, 6)

# vertex labels as strings
labels = string.(1:nv(g))

# hide vertices with even indices
vis = [isodd(i) for i in 1:nv(g)]
node_color = [vis[i] ? :black : RGBAf(0, 0, 0, 0) for i in 1:nv(g)]
node_size = [vis[i] ? 50 : 0 for i in 1:nv(g)]
label_color = [vis[i] ? :white : RGBAf(0, 0, 0, 0) for i in 1:nv(g)]

fig = Figure(size = (900, 500))
ax = Axis(
    fig[1, 1];
    xgridvisible = false,
    ygridvisible = false,
    xticksvisible = false,
    yticksvisible = false,
    xticklabelsvisible = false,
    yticklabelsvisible = false,
)
hidespines!(ax)

graphplot!(
    ax,
    g;
    layout = Spring(C = 50.0),
    nlabels = labels,
    nlabels_color = label_color,
    nlabels_size = 18,
    nlabels_offset = Point2f(0, 0),
    nlabels_align = (:center, :center),
    node_color = node_color,
    node_size = node_size,
    arrow_shift = :end,
    arrow_size = 12,
    edge_width = 2,
    edge_curvature = 0.9,
)

display(fig)

using GLMakie
using GraphMakie
using Graphs
g = wheel_graph(10)
f, ax, p = graphplot(g, edge_width = [3 for i in 1:ne(g)], node_size = [10 for i in 1:nv(g)])

deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))
register_interaction!(ax, :edrag, EdgeDrag(p))

deregister_interaction!(ax, :nhover)
deregister_interaction!(ax, :ehover)
deregister_interaction!(ax, :ndrag)
deregister_interaction!(ax, :edrag)


using Makie.Colors

function action(idx, event, axis)
    println("Clicked on node index: $idx")
    p.node_color[][idx] = rand(RGB)
    p.node_color[] = p.node_color[]
    p.node_size[][idx] = p.node_size[][idx] + 5
    p.node_size[] = p.node_size[]
    println("New color: $(p.node_color[][idx]), New size: $(p.node_size[][idx])")
end

g = wheel_digraph(10)
f, ax, p = graphplot(g, node_size = [30 for i in 1:nv(g)], node_color = [colorant"red" for i in 1:nv(g)])

deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :nodeclick, NodeClickHandler(action))