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

const stop_flag = Base.Threads.Atomic{Bool}(false)

process_strong_c3_combinations_by_combination_size(
    GCATCodes.ALL_CODONS,
    1,
    "files/results/test1.txt",
    "files/checkpoints/test1_cp.txt",
    stop_flag;
    show_debug = true,
)


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

_get_processed_count_from_combination([1, 2, 3])

for i in 7:10
    println("--------------------------------- Combination size: $i ---------------------------------")
    checkpoint = _load_strong_c3_checkpoint("files/checkpoints/test$(i)_cp.txt")
    comb = checkpoint.current_combination
    if length(comb) > i
        # set comb as last combination of size i
        comb = collect((60 - i + 1):60)
        println("SET comb to last combination for size $i: $comb")
    end
    res = _get_processed_count_from_combination(comb)
    # get amount of lines in result file
    lines = readlines("files/results/test$(i).txt")
    processed_count = res
    strong_c3_count = length(lines)
    not_strong_c3_count = processed_count - strong_c3_count
    println(
        "processed_count: $processed_count, strong_c3_count: $strong_c3_count, not_strong_c3_count: $not_strong_c3_count",
    )
end

process_strong_c3_combinations_by_combination_size(
    GCATCodes.ALL_CODONS,
    1,
    "files/test.txt",
    "files/test_cp.txt",
    stop_flag;
    show_debug = true,
)


lengths = 1:20
tasks = [
    @spawn begin
        k = len
        res = "files/tests/strong_c3_cs$(k).txt"
        chk = "files/tests/strong_c3_cs$(k)_cp.txt"
        process_strong_c3_combinations_by_combination_size(
            GCATCodes.ALL_CODONS,
            k,
            res,
            chk,
            false,               # kein Checkpoint laden
            Atomic{Bool}(false); # eigener Cancel-Schalter
            show_debug = false,
        )
    end for len in lengths
]

fetch.(tasks)  # wartet, bis alle fertig sind


using Base.Threads: Atomic
using Base.Threads
stopflag = Atomic{Bool}(false)

process_strong_c3_combinations_by_combination_size(
    GCATCodes.ALL_CODONS,
    6,
    "files/tests/strong_c3_cs5.txt",
    "files/tests/strong_c3_cs5_cp.txt",
    false,
    stopflag;
    show_debug = false,
)

stopflag[] = true

const MAX_LENGTH::Int = 20
const STRONG_C3_RESULTS_PATH = "files/results/strong_c3_codon_combinations.txt"
const STRONG_C3_CHECKPOINT_PATH = "files/checkpoints/strong_c3_checkpoint.txt"
const TEST_RESULTS_PATH = "files/results/test_strong_c3_codon_combinations.txt"
const TEST_CHECKPOINT_PATH = "files/checkpoints/test_strong_c3_checkpoint.txt"

# TEST
# start fresh
const TEST_LENGTH = 4
process_strong_c3_combinations(GCATCodes.ALL_CODONS, 2, TEST_RESULTS_PATH, TEST_CHECKPOINT_PATH, false)

# resume
process_strong_c3_combinations(
    GCATCodes.ALL_CODONS,
    TEST_LENGTH,
    TEST_RESULTS_PATH,
    TEST_CHECKPOINT_PATH,
    true,
)

@benchmark process_strong_c3_combinations(
    GCATCodes.ALL_CODONS,
    TEST_LENGTH,
    TEST_RESULTS_PATH,
    TEST_CHECKPOINT_PATH,
    false,
) samples = 5


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

# -------------------------------------------------- TESTING PARALLELISM --------------------------------------------------
using BenchmarkTools
using Base.Threads: Atomic
using Base.Threads

Threads.nthreads()


function test1(
    codons::Vector{LongDNA{4}},
    max_combination_size_length::Int,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    load_checkpoint::Bool;
    show_debug::Bool = false,
    cancel::Atomic{Bool} = Atomic{Bool}(false),
)
    # load from checkpoint or start from scratch
    check_point = load_checkpoint ? _load_strong_c3_checkpoint(checkpoint_path) : nothing
    current_combination = collect(1:max_combination_size_length)
    processed_count = load_checkpoint ? check_point.processed_count : 0
    strong_c3_count = load_checkpoint ? check_point.strong_c3_count : 0
    not_strong_c3_count = load_checkpoint ? check_point.not_strong_c3_count : 0

    # adjust variables based on load_checkpoint
    write_mode = load_checkpoint ? "a" : "w"

    # get length of codon set
    length_codon_set = length(codons)
    # disable default SIGINT handler to allow custom handling
    try
        open(results_path, write_mode) do result_out
            # iterate over all combination sizes
            for combination_size in max_combination_size_length:max_combination_size_length
                while true
                    # allow interrupting the process
                    cancel[] && throw(InterruptException())

                    codon_set = codons[current_combination]

                    # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
                    if _contains_codon_rotation(codon_set)
                        not_strong_c3_count += 1
                    else
                        # check strong C3
                        data = CodonGraphData(codon_set)
                        construct_graph_data!(data; show_debug = false)
                        if is_strong_c3(data; show_debug = false)
                            codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                            println(
                                result_out,
                                "Strong C3: $codon_combination_string with size $(length(codon_set))",
                            )
                            strong_c3_count += 1
                        else
                            not_strong_c3_count += 1
                        end
                    end

                    processed_count += 1

                    # get next combination or break if none left
                    if !_increment_codon_set_combination!(current_combination, length_codon_set)
                        current_combination = collect(1:(combination_size + 1))
                        break
                    end
                end
            end
        end
    catch err
        err isa InterruptException || rethrow(err)
    finally
    end
end

function test2(
    codons::Vector{LongDNA{4}},
    max_combination_size_length::Int,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    load_checkpoint::Bool;
    show_debug::Bool = false,
    cancel::Atomic{Bool} = Atomic{Bool}(false),
)
    # load from checkpoint or start from scratch
    check_point = load_checkpoint ? _load_strong_c3_checkpoint(checkpoint_path) : nothing
    current_combination = collect(1:max_combination_size_length)
    processed_count = load_checkpoint ? check_point.processed_count : 0
    strong_c3_count = load_checkpoint ? check_point.strong_c3_count : 0
    not_strong_c3_count = load_checkpoint ? check_point.not_strong_c3_count : 0

    # adjust variables based on load_checkpoint
    write_mode = load_checkpoint ? "a" : "w"

    # get length of codon set
    length_codon_set = length(codons)
    # disable default SIGINT handler to allow custom handling
    try
        open(results_path, write_mode) do result_out
            # iterate over all combination sizes
            for combination_size in max_combination_size_length:max_combination_size_length
                while true
                    # allow interrupting the process
                    cancel[] && throw(InterruptException())

                    codon_set = codons[current_combination]

                    # check strong C3
                    data = CodonGraphData(codon_set)
                    construct_graph_data!(data; show_debug = false)
                    if is_strong_c3(data; show_debug = false)
                        codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                        println(
                            result_out,
                            "Strong C3: $codon_combination_string with size $(length(codon_set))",
                        )
                        strong_c3_count += 1
                    else
                        not_strong_c3_count += 1
                    end

                    processed_count += 1

                    # get next combination or break if none left
                    if !_increment_codon_set_combination!(current_combination, length_codon_set)
                        current_combination = collect(1:(combination_size + 1))
                        break
                    end
                end
            end
        end
    catch err
        err isa InterruptException || rethrow(err)
    finally
    end
end


stopflag = Atomic{Bool}(false)

test1(
    GCATCodes.ALL_CODONS,
    3,
    "files/tests/test1_results.txt",
    "files/tests/test1_cp.txt",
    false;
    cancel = stopflag,
)

test2(
    GCATCodes.ALL_CODONS,
    3,
    "files/tests/test2_results.txt",
    "files/tests/test2_cp.txt",
    false;
    cancel = stopflag,
)


TEST_LENGTH2 = 4
# @benchmark test1(
#     $GCATCodes.ALL_CODONS,
#     $TEST_LENGTH2,
#     "files/tests/test1_results.txt",
#     "files/tests/test1_cp.txt",
#     false;
#     cancel = stopflag,
# ) samples = 5
# @benchmark test2(
#     $GCATCodes.ALL_CODONS,
#     $TEST_LENGTH2,
#     "files/tests/test2_results.txt",
#     "files/tests/test2_cp.txt",
#     false;
#     cancel = stopflag,
# ) samples = 5


t1 = Threads.@spawn begin
    @benchmark test1(
        $GCATCodes.ALL_CODONS,
        $TEST_LENGTH2,
        "files/tests/test1_results.txt",
        "files/tests/test1_cp.txt",
        false;
        cancel = stopflag,
    ) samples = 5
end

t2 = Threads.@spawn begin
    @benchmark test2(
        $GCATCodes.ALL_CODONS,
        $TEST_LENGTH2,
        "files/tests/test2_results.txt",
        "files/tests/test2_cp.txt",
        false;
        cancel = stopflag,
    ) samples = 5
end

t3 = Threads.@spawn begin
    @benchmark test1(
        $GCATCodes.ALL_CODONS,
        8,
        "files/tests/test3_results.txt",
        "files/tests/test3_cp.txt",
        false;
        cancel = stopflag,
    ) samples = 5
end

t4 = Threads.@spawn begin
    @benchmark test2(
        $GCATCodes.ALL_CODONS,
        8,
        "files/tests/test4_results.txt",
        "files/tests/test4_cp.txt",
        false;
        cancel = stopflag,
    ) samples = 5
end

fetch.((t1, t2))
fetch.((t3, t4))

for (i, t) in enumerate(tasks)
    @info "task" idx = i done = istaskdone(t) failed = istaskfailed(t) state = t.state
end

# compare if 2 files are identical
same = read("files/tests/test1_results.txt") == read("files/tests/test2_results.txt")
same = read("files/tests/test3_results.txt") == read("files/tests/test4_results.txt")



# abbrechen
stopflag[] = true
foreach(t -> istaskdone(t) || istaskfailed(t) || schedule(t, InterruptException()), tasks)
# warten und Interrupts schlucken
for t in tasks
    try
        fetch(t)
    catch err
        err isa InterruptException || rethrow(err)
    end
end




function test3(
    codons::Vector{LongDNA{4}},
    max_combination_size_length::Int,
    results_path::AbstractString,
    checkpoint_path::AbstractString,
    load_checkpoint::Bool;
    show_debug::Bool = false,
    cancel::Atomic{Bool} = Atomic{Bool}(false),
)
    # load from checkpoint or start from scratch
    check_point = load_checkpoint ? _load_strong_c3_checkpoint(checkpoint_path) : nothing
    current_combination = collect(1:max_combination_size_length)
    processed_count = load_checkpoint ? check_point.processed_count : 0
    strong_c3_count = load_checkpoint ? check_point.strong_c3_count : 0
    not_strong_c3_count = load_checkpoint ? check_point.not_strong_c3_count : 0

    # adjust variables based on load_checkpoint
    write_mode = load_checkpoint ? "a" : "w"

    # get length of codon set
    length_codon_set = length(codons)
    # disable default SIGINT handler to allow custom handling
    try
        open(results_path, write_mode) do result_out
            # iterate over all combination sizes
            for combination_size in max_combination_size_length:max_combination_size_length
                while true
                    # allow interrupting the process
                    cancel[] && throw(InterruptException())

                    codon_set = codons[current_combination]

                    # skip combinations that contain N1N2N3, N2N3N1 and N3N1N2 for some codon
                    if _contains_codon_rotation(codon_set)
                        not_strong_c3_count += 1
                    else
                        # check strong C3
                        data = CodonGraphData(codon_set)
                        construct_graph_data!(data; show_debug = false)
                        if is_strong_c3(data; show_debug = false)
                            codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
                            println(
                                result_out,
                                "Strong C3: $codon_combination_string with size $(length(codon_set))",
                            )
                            strong_c3_count += 1
                        else
                            not_strong_c3_count += 1
                        end
                    end

                    processed_count += 1

                    # get next combination or break if none left
                    if !_increment_codon_set_combination!(current_combination, length_codon_set)
                        current_combination = collect(1:(combination_size + 1))
                        break
                    end
                end
            end
        end
    catch err
        err isa InterruptException || rethrow(err)
    finally
    end
end

min_length3 = 1
max_length3 = 10

tasks = [
    @spawn begin
        tid = Threads.threadid()
        println("Thread $tid started job $i")
        try
            test3(
                GCATCodes.ALL_CODONS,
                i,
                "files/tests/test_results_task_$i.txt",
                "files/tests/test_checkpoint_task_$i.txt",
                false;
                cancel = stopflag,
            )
            println("Thread $tid finished job $i")
        catch err
            if err isa InterruptException
                println("Thread $tid interrupted job $i")
            else
                rethrow(err)
            end
        end
    end for i in min_length3:max_length3
]


# abbrechen
stopflag[] = true
foreach(t -> istaskdone(t) || istaskfailed(t) || schedule(t, InterruptException()), tasks)
# warten und Interrupts schlucken
for t in tasks
    try
        fetch(t)
    catch err
        err isa InterruptException || rethrow(err)
    end
end

Threads.@threads for _ in 1:Threads.nthreads()
    @info "alive" thread = Threads.threadid()
end

for (i, t) in enumerate(tasks)
    @info "task" idx = i state = t.state failed = istaskfailed(t) done = istaskdone(t)
end


codon_set = LongDNA{4}.(GCATCodes.ALL_CODONS[5:15])
count = 0


using Random

num = 2
for i in num:num
    idxs = randperm(length(GCATCodes.ALL_CODONS))[1:i]
    codon_set = LongDNA{4}.(GCATCodes.ALL_CODONS[idxs])
    println("i: $i")
    t1 = @benchmark ab1(count, codon_set)
    t2 = @benchmark ab2(count, codon_set)
    m1 = median(t1).time
    m2 = median(t2).time

    if m1 > m2
        println("ab2 faster: m2: $(m2)ns, m1: $(m1)ns, ratio: $(round(m1 / m2, digits=2))")
    else
        println("ab1 faster: m2: $(m2)ns, m1: $(m1)ns, ratio: $(round(m2 / m1, digits=2))")
    end
end

function ab1(count, codon_set)
    open("files/tests/ab1_results.txt", "w") do result_out
        data = CodonGraphData(codon_set)
        construct_graph_data!(data)
        if is_strong_c3(data; show_debug = false)
            codon_combination_string = join("\"" .* string.(codon_set) .* "\"", ", ")
            println(result_out, "Strong C3: $codon_combination_string with size $(length(codon_set))")
            count += 1
        else
            count += 1
        end
    end
end

function ab2(count, codon_set)
    if _contains_codon_rotation(codon_set)
        count += 1
    end
end

ab1(count, codon_set)
ab2(count, codon_set)
