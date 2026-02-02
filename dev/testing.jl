using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs
using BioSequences
using NetworkLayout
using Base.Threads
using BenchmarkTools

# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
global_logger(ConsoleLogger(Logging.Info)) # deactivate

open("files/results/test.txt", "r") do res
    counter = 0
    for line in eachline(res)
        if counter == 100_000
            break
        end
        codon_set = result_to_codon_set(line)
        # println("Codon set: $codon_set")
        counter += 1
    end
end

open("files/results/test.txt", "w") do in
    open("files/results/result_8.txt", "r") do out
        for line in eachline(out)
            println(in, line)
        end
    end
end


open("files/results/result_5.txt", "r") do res
    counter = 0
    false_positives = 0
    for line in eachline(res)
        # if counter == 10
        #     break
        # end
        counter += 1
        if counter % 100000 == 0
            println("counter: $counter")
        end

        codon_set = result_to_codon_set(line)
        data = CodonGraphData(codon_set)
        construct_graph_data!(data)
        _expand_graph(data)

        if get_max_cycle_length(data.graph; show_debug = false) > 2
            println("Max cycle length > 2 found!")
            println("codon set: $codon_set in line $counter is FALSE POSITIVE")
            false_positives += 1
        end
    end
    println("FINISHED, counter: $counter, false_positives: $false_positives")
end

function subsets(codon_set::Vector{LongDNA{4}})
    n = length(codon_set)
    result = Vector{Vector{LongDNA{4}}}()
    for i in 1:(2^n - 1)
        subset = Vector{LongDNA{4}}()
        for j in 1:n
            if (i >> (j - 1)) & 1 == 1
                push!(subset, codon_set[j])
            end
        end
        push!(result, subset)
    end
    return result
end


set = LongDNA{4}.(["ACG", "GCC", "GGC"])
data = CodonGraphData(set)
construct_graph_data!(data)
show_codon_graph(data)







const stop_flag = Base.Threads.Atomic{Bool}(false)
codons = GCATCodes.ALL_CODONS


function ab1()
    s = 0
    for i in 1:iters
        if cancel[] && (s += 1)
        end
    end
    s
end

function ab2(c)
    s = 0
    for i in 1:iters
        c -= 1
        if c == 0
            cancel[] && (s += 1)
            c = 100_000
        end
    end
    s
end

cancel = Base.Threads.Atomic{Bool}(false)
iters = 100_000_000
c = 100_000

ab1()
ab2(c)

@btime ab1() samples = 500 evals = 1
@btime ab2(c) samples = 500 evals = 1


# check if there is an empty line in a file
function check_empty_lines_in_file(file_path::String)
    open(file_path, "r") do f
        line_number = 0
        for line in eachline(f)
            line_number += 1
            if isempty(strip(line))
                println("Empty line found at line number: $line_number")
            end
        end
    end
end

check_empty_lines_in_file("files/results/result_7.txt")
check_empty_lines_in_file("files/results/result_8.txt")
check_empty_lines_in_file("files/results/result_9.txt")
check_empty_lines_in_file("files/results/result_10.txt")
check_empty_lines_in_file("files/results/result_11.txt")
check_empty_lines_in_file("files/results/result_12.txt")
check_empty_lines_in_file("files/results/result_13.txt")
check_empty_lines_in_file("files/results/result_14.txt")
check_empty_lines_in_file("files/results/result_15.txt")
check_empty_lines_in_file("files/results/result_16.txt")
check_empty_lines_in_file("files/results/result_17.txt")
check_empty_lines_in_file("files/results/result_18.txt")
check_empty_lines_in_file("files/results/result_19.txt")
check_empty_lines_in_file("files/results/result_20.txt")

# function to get combination from codon set
function _get_combination_from_codon_set(codon_set::Vector{LongDNA{4}}, codons::Vector{LongDNA{4}})
    combination = Vector{Int}()
    codon_set_set = Set(codon_set)
    @inbounds for (i, c) in enumerate(codons)
        if c in codon_set_set
            push!(combination, i)
        end
    end
    return combination
end

codon_set = LongDNA{4}.(["AGC", "ATA", "ATG", "CAC", "CTC", "GTG", "TTG"])
codon_set2 = LongDNA{4}.(["AGC", "ATA", "ATG", "CAC", "CTC", "TAC", "TAG"])
println(_get_combination_from_codon_set(codon_set, codons))
println(_get_combination_from_codon_set(codon_set2, codons))

combination = _get_combination_from_codon_set(codon_set, codons)
combination2 = _get_combination_from_codon_set(codon_set2, codons)

for _ in 1:1
    counter = 0
    _increment_codon_set_combination!(combination, 60)
    while combination != combination2
        counter += 1
        codon_set = codons[combination]
        println("combination: $combination")
        println("codon_set: $codon_set")
        data = CodonGraphData(codon_set)
        construct_graph_data!(data; show_debug = false)
        if is_strong_c3(data; show_debug = false)
            println("strong C3: $codon_set")
        end
        _increment_codon_set_combination!(combination, 60)
    end
    println(counter)
end








# function to get processed count from current combination
function _get_processed_count_from_combination(comb::Vector{Int}; n::Int = 60)
    k = length(comb)
    rank = 0
    prev = 0
    @inbounds for (i, ci) in enumerate(comb)
        for x in (prev + 1):(ci - 1)
            rank += binomial(n - x, k - i)
        end
        prev = ci
    end
    return rank
end
# ---------------------------------------------- FUNCTIONS ----------------------------------------------
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

const MICRO_ITERS = 1_000_000
codons = GCATCodes.ALL_CODONS
rotation_masks = GCATCodes._build_rotation_masks(codons)

index_map = Dict{LongDNA{4}, Int}()
@inbounds for (i, c) in enumerate(codons)
    index_map[c] = i
end

base_idx = 1
base = codons[base_idx]
idx1 = index_map[left_shift_codon(base, 1)]
idx2 = index_map[left_shift_codon(base, 2)]
extra_idx = base_idx == 1 ? 2 : 1
if extra_idx == idx1 || extra_idx == idx2
    extra_idx = 3
end
combination = [base_idx, idx1, idx2, extra_idx]

function _reject_plain_once(codons, combination)
    codon_set = codons[combination]
    return _contains_codon_rotation(codon_set)
end

function _reject_mask_once(rotation_masks, combination)
    comb_mask = GCATCodes._combination_to_mask(combination)
    return GCATCodes._mask_contains_codon_rotation(comb_mask, combination, rotation_masks)
end

function _micro_reject_plain(codons, combination, iters)
    hits = 0
    for _ in 1:iters
        _reject_plain_once(codons, combination) && (hits += 1)
    end
    return hits
end

function _micro_reject_mask(rotation_masks, combination, iters)
    hits = 0
    for _ in 1:iters
        _reject_mask_once(rotation_masks, combination) && (hits += 1)
    end
    return hits
end

trial_plain = @benchmark _micro_reject_plain($codons, $combination, $MICRO_ITERS) samples = 500 evals = 1
trial_mask =
    @benchmark _micro_reject_mask($rotation_masks, $combination, $MICRO_ITERS) samples = 500 evals = 1
m_plain = median(trial_plain).time
m_mask = median(trial_mask).time
println("plain: $(m_plain) ns for $MICRO_ITERS iterations")
println("mask : $(m_mask) ns for $MICRO_ITERS iterations")
if m_mask < m_plain
    println("speedup: $(round(m_plain / m_mask, digits = 2))x (mask faster)")
else
    println("speedup: $(round(m_mask / m_plain, digits = 2))x (plain faster)")
end



# inkrementelle strong‑C3-Suche: wächst nur von bereits gefundenen Treffern weiter
function run_incremental_strong_c3!(
    codons::Vector{LongDNA{4}} = GCATCodes.ALL_CODONS;
    max_size::Int = 20,
    results_dir::AbstractString = "files/results",
    show_debug::Bool = false,
)
    rot_masks = GCATCodes._get_rotation_masks(codons)
    frontier = Vector{Vector{Int}}()

    # Größe 1 seeds
    for i in 1:length(codons)
        combo = [i]
        _combo_is_strong_c3!(combo, codons, rot_masks, results_dir, 1; show_debug) && push!(frontier, combo)
    end

    # Größen 2..max_size
    for size in 2:max_size
        isempty(frontier) && break
        next_frontier = Vector{Vector{Int}}()

        for base in frontier
            last = base[end]
            for j in (last + 1):length(codons)           # aufsteigend -> keine Duplikate
                combo = [base...; j]
                mask = GCATCodes._combination_to_mask(combo)
                GCATCodes._mask_contains_rotation(combo, mask, rot_masks) && continue
                _combo_is_strong_c3!(combo, codons, rot_masks, results_dir, size; show_debug) &&
                    push!(next_frontier, combo)
            end
        end

        frontier = next_frontier
    end
end

# Einzel-Kombi prüfen und ggf. in Ergebnisdatei schreiben
function _combo_is_strong_c3!(
    combo::Vector{Int},
    codons::Vector{LongDNA{4}},
    rot_masks::Vector{UInt64},
    results_dir::AbstractString,
    size::Int;
    show_debug::Bool = false,
)
    data = CodonGraphData(codons[combo])
    construct_graph_data!(data; show_debug = false)
    if is_strong_c3(data; show_debug = show_debug)
        path = joinpath(results_dir, "result_inc_$size.txt")
        open(path, "a") do io
            codon_str = join("\"" .* string.(codons[combo]) .* "\"", ", ")
            println(io, "Strong C3: $codon_str with size $size")
        end
        return true
    end
    return false
end
