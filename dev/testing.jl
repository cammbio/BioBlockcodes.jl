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
# function to add N₂ and N₃N₁ vertices and connect edges for whole codon set
function add_n2_n3n1_vertices_and_edges!(data::CodonGraphData; show_debug::Bool = false)
    for codon in data.codon_set
        n2 = string(codon[2])
        n3n1 = string(codon[3], codon[1])
        add_vertex_by_label!(data, n2, show_debug = false)
        add_vertex_by_label!(data, n3n1, show_debug = false)
        add_edge_by_label!(data, n2, n3n1, show_debug = false)
        add_edge_by_label!(data, n3n1, n2, show_debug = false)
    end
end

# function to add N₂ and N₃N₁ vertices and connect edges for a single codon
function add_n2_n3n1_vertices_and_edges_for_codon!(
    data::CodonGraphData,
    codon::LongDNA{4};
    show_debug::Bool = false,
)
    n2 = string(codon[2])
    n3n1 = string(codon[3], codon[1])
    add_vertex_by_label!(data, n2, show_debug = false)
    add_vertex_by_label!(data, n3n1, show_debug = false)
    add_edge_by_label!(data, n2, n3n1, show_debug = false)
    add_edge_by_label!(data, n3n1, n2, show_debug = false)
end

function check_properties(data::CodonGraphData)
    println("Checking properties of codon graph with codon set: $(data.codon_set)")
    println("is_c3: $(is_c3(data; show_debug = false))")
    println("is_circular: $(is_circular(data.graph; show_debug = false))")
    println("is_comma_free: $(is_comma_free(data.graph; show_debug = false))")
    println("is_self_complementary: $(is_self_complementary(data; show_debug = false))")
end


# function to turn file line to codon set
function line_to_codon_set(line::String)
    codon_array = split(line, ", ")
    codon_array = replace.(codon_array, "\"" => "")
    codon_set = LongDNA{4}.(codon_array)
    return codon_set
end


# generate all combinations of a codon set by a specific size
function all_combinations_per_size(codon_set::Vector{LongDNA{4}}, combination_size::Int)
    length_codon_set = length(codon_set)
    combination_size < 0 && throw(ArgumentError("combination_size must be ≥ 0"))
    combination_size == 0 && return [LongDNA{4}[]]
    combination_size > length_codon_set && return Vector{Vector{LongDNA{4}}}()

    combos = Vector{Vector{LongDNA{4}}}()
    combination = collect(1:combination_size)

    push!(combos, codon_set[combination])  # erste Kombination

    while _next_combination!(combination, combination_size, length_codon_set)
        push!(combos, codon_set[combination])
    end

    return combos
end


# 
function _next_combination!(combination::Vector{LongDNA{4}}, combination_size::Int, length_codon_set::Int)
    for i in combination_size:-1:1
        if combination[i] != i + length_codon_set - combination_size
            combination[i] += 1
            for j in (i + 1):combination_size
                combination[j] = combination[j - 1] + 1
            end
            return true
        end
    end
    return false
end



# function to get each possible combination of codons from a codon set
function get_codon_combinations(codon_set::Vector{LongDNA{4}})
    n = length(codon_set)
    codon_combinations = Vector{Vector{LongDNA{4}}}()
    for i in 1:n
        combs = collect(combinations(codon_set, i))
        append!(codon_combinations, combs)
    end
    return codon_combinations
end


# -------------------------------------------------- TESTING --------------------------------------------------
codon_set = LongDNA{4}.(["AAC", "GTT", "AAG", "CTT", "AAT", "ATT", "ACC"])
data_original = CodonGraphData(codon_set; plot_title = "Original")
construct_graph_data!(data_original; show_debug = false)

data_expanded = CodonGraphData(codon_set; plot_title = "Expanded")
construct_graph_data!(data_expanded; show_debug = false)
add_n2_n3n1_vertices_and_edges!(data_expanded; show_debug = false)

a = get_cycles_difference(data_expanded, data_original; show_debug = false)
println(a)


# read file per line with all 216 maximal self-complementary C3 codon_set and write a new file where each line is written as Array
in_path = "files/216_maximal_self_complementary_c3_codes.txt"
maximal_c3_codons_path = "files/216_maximal_self_complementary_c3_codes_array.txt"
open(maximal_c3_codons_path, "w") do out
    for line in eachline(in_path)
        codon_array = split(line)
        formatted = join("\"" .* codon_array .* "\"", ", ")
        println(out, formatted)
    end
end



open(maximal_c3_codons_path, "r") do f
    codon_set_size_count = Dict{Int, Int}()
    for line in eachline(f)
        codon_set = line_to_codon_set(line)
        codon_set_growing = LongDNA{4}[]
        for codon in codon_set
            push!(codon_set_growing, codon)
            # create codon graph data
            data = CodonGraphData(codon_set_growing; plot_title = "Codon set: $codon_set_growing")
            construct_graph_data!(data; show_debug = false)
            # manually add N₂ and N₃N₁ vertices and connect edges
            add_n2_n3n1_vertices_and_edges!(data; show_debug = false)
            # show_graph(data; show_debug = false)

            # check if codon graph is C3 and comma-free
            if is_comma_free(data.graph; show_debug = false) && is_c3(data; show_debug = false)
                # show_graph(data; show_debug = false)
                codon_set_length = length(codon_set_growing)
                codon_set_size_count[codon_set_length] = get(codon_set_size_count, codon_set_length, 0) + 1
                open("files/c3_comma_free_codes.txt", "a") do out
                    redirect_stdout(out) do
                        codon_set_growing_string = join("\"" .* string.(codon_set_growing) .* "\"", ", ")
                        println("""Current codon set $codon_set_growing_string is C3 and comma-free.
                        codon set length: $(length(codon_set_growing))""")
                        println()
                    end
                end
            end
        end
    end
    open("files/codon_set_size_count_after_adding_vertices_and_edges.txt", "w") do out
        redirect_stdout(out) do
            println("Size count for C3 and comma-free codon sets after adding vertices and edges:")
            for (size, count) in sort(collect(codon_set_size_count))
                println("Size $size: $count")
            end
        end
    end

    println("Growing codon set analysis complete. Results written to files/c3_comma_free_codes.txt")
end


# group codon sets by maximal cycle length
open(maximal_c3_codons_path, "r") do f
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
open(maximal_c3_codons_path, "r") do f
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
open(maximal_c3_codons_path, "r") do f
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


# check for each codon in each codon set if adding N₂ and N₃N₁ vertices and edges (per codon) keeps the graph C3 
open(maximal_c3_codons_path, "r") do f
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

                    add_n2_n3n1_vertices_and_edges_for_codon!(data, codon; show_debug = false)
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


# check for each codon set which combination of codons keeps the graph C3 after adding N₂ and N₃N₁ vertices and edges
open(maximal_c3_codons_path, "r") do f
    open("files/codon_combinations_keeping_c3_property.txt", "w") do out
        redirect_stdout(out) do
            for line in eachline(f)
                # construct graph data
                codon_set = line_to_codon_set(line)
                println("Codon set: $line")
                n = length(codon_set)
                # check all combinations of codons
                for i in 1:n
                    codon_combinations = collect(combinations(codon_set, i))
                    for codon_combination in codon_combinations
                        data = CodonGraphData(codon_set; plot_title = "Codon set: $codon_set")
                        construct_graph_data!(data; show_debug = false)

                        for codon in codon_combination
                            add_n2_n3n1_vertices_and_edges_for_codon!(data, codon; show_debug = false)
                        end

                        result = is_c3(data; show_debug = false)
                        if result
                            codon_combination_string = join("\"" .* string.(codon_combination) .* "\"", ", ")
                            println("Codon combination keeping C3 property: $codon_combination_string")
                        end
                    end
                end
            end
        end
    end
    println(
        "Codon combination analysis complete. Results written to files/codon_combinations_keeping_c3_property.txt",
    )
end


codon_set =
    LongDNA{4}.(["CAA", "TTG", "CAC", "GTG", "CAG", "CTG", "CTC", "GAG", "GAA", "TTC", "GAC", "GTC", "GCC"])

println(all_combinations_per_size(codon_set, 12))
println(length(all_combinations_per_size(codon_set, 12)))


codon_set =
    LongDNA{4}.(["CAA", "TTG", "CAC", "GTG", "CAG", "CTG", "CTC", "GAG", "GAA", "TTC", "GAC", "GTC", "GCC"])
data = CodonGraphData(codon_set)
construct_graph_data!(data; show_debug = false)
show_graph(data; show_debug = false)
add_n2_n3n1_vertices_and_edges!(data; show_debug = false)
cycles = get_all_cycles(data; show_debug = false)
println(typeof(cycles))
for cycle in cycles
    println(cycle)
end


# compare α₁ and α₂ to original graph
data_original = CodonGraphData(codon_set)
construct_graph_data!(data_original; show_debug = false)
data_alpha1 = CodonGraphData(left_shift_codon_set(codon_set, 1))
construct_graph_data!(data_alpha1; show_debug = false)
data_alpha2 = CodonGraphData(left_shift_codon_set(codon_set, 2))
construct_graph_data!(data_alpha2; show_debug = false)
show_graph(data_original; show_debug = false)
show_graph(data_alpha1; show_debug = false)
show_graph(data_alpha2; show_debug = false)

# check which vertices are common between original and α₁ and α₂ graphs
common_vertices_alpha1 = intersect(keys(data_original.vertex_index), keys(data_alpha1.vertex_index))
common_vertices_alpha2 = intersect(keys(data_original.vertex_index), keys(data_alpha2.vertex_index))
println("Common vertices between original and α₁ graph: $common_vertices_alpha1")
println("Common vertices between original and α₂ graph: $common_vertices_alpha2")
println("data_original.vertex_index: $(keys(data_original.vertex_index))")
println("data_alpha1.vertex_index: $(keys(data_alpha1.vertex_index))")
println("data_alpha2.vertex_index: $(keys(data_alpha2.vertex_index))")

# function to compare two codon graph datas
function compare_codon_graph_data(data1::CodonGraphData, data2::CodonGraphData)
    vertices1 = keys(data1.vertex_index)
    vertices2 = keys(data2.vertex_index)
    common_vertices = intersect(vertices1, vertices2)
    println("Common vertices: $common_vertices")

    edges1 = collect(edges(data1.graph))
    edges2 = collect(edges(data2.graph))
    common_edges = intersect(edges1, edges2)
    println("Common edges: $common_edges")
end

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

open("files/test123", "r") do f
    for line in eachline(f)
        println(line)
        expr = Meta.parse(line)
        codon_set = eval(expr)
        data = CodonGraphData(codon_set; plot_title = "Codon set: $line")
        construct_graph_data!(data; show_debug = false)
        show_graph(data; show_debug = false)
    end
end














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