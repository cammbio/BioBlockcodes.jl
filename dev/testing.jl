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


# read file per line with all 216 maximal self-complementary C3 codon_set and write a new file where each line is written as Array
in_path = "files/216_maximal_self_complementary_c3_codes.txt"
out_path = "files/216_maximal_self_complementary_c3_codes_array.txt"
open(out_path, "w") do out
    for line in eachline(in_path)
        codon_array = split(line)
        formatted = join("\"" .* codon_array .* "\"", ", ")
        println(out, formatted)
    end
end

open(out_path, "r") do f
    # write all println output to file
    println("Running loop over all codon sets in file and performing analysis...")
    open("files/demo_output.txt", "w") do out
        redirect_stdout(out) do

            counter = 0
            for line in eachline(f)
                # break for debugging
                counter += 0 # set to 0 to run all
                if counter > 1
                    break
                end

                # turn line into codon set
                codon_set = line_to_codon_set(line)
                # create codon graph data
                data_original = CodonGraphData(codon_set; plot_title = "Codon set: $codon_array before")
                construct_graph_data!(data_original; show_debug = false)


                # create second CodonGraphData from original codon_set to compare cycles
                data_adjusted = CodonGraphData(codon_set; plot_title = "Codon set: $codon_array after")
                construct_graph_data!(data_adjusted; show_debug = false)
                # manually add N₂ and N₃N₁ vertices and connect edges
                add_n2_n3n1_vertices_and_edges!(data_adjusted, codon_set; show_debug = false)

                # only proceed if adjusted graph is still C3
                if is_c3(data_adjusted; show_debug = false)

                    # display initial graph
                    show_graph(data_original; show_debug = false)
                    # check mathematical properties of graph before adding new vertices/edges
                    check_properties(data_original)
                    # display all cycles before adding new vertices/edges
                    display_all_cycles(data_original; show_debug = false)
                    # get duplicate cycles before adding new vertices/edges
                    display_duplicate_cycles(data_original; show_debug = false)
                    # get cycles of all lengths before adding new vertices/edges
                    check_cycles(data_original)


                    # display updated graph
                    show_graph(data_adjusted; show_debug = false)
                    # check mathematical properties of graph after adding new vertices/edges
                    check_properties(data_adjusted)
                    # display all cycles after adding new vertices/edges
                    # println_to_file("files/demo_output.txt", () -> display_all_cycles(data_adjusted; show_debug = false))
                    display_all_cycles(data; show_debug = false)
                    # get duplicate cycles after adding new vertices/edges
                    display_duplicate_cycles(data_adjusted; show_debug = false)
                    # get cycles of all lengths after adding new vertices/edges
                    check_cycles(data_adjusted)
                    # get cycles with length 2
                    get_cycles_of_length(data_adjusted, 2; show_debug = false)


                    # get difference in cycles between original and adjusted graphs
                    get_cycles_difference(data_original, data_adjusted; show_debug = false)
                end
            end
        end
    end
    println("Analysis complete. Results written to files/demo_output.txt")
end


open(out_path, "r") do f
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
                        # println()
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
open(out_path, "r") do f
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
        println(
            "Maximal cycle length: $maximal_cycle_length -----------------------------------------------------------------------------------------",
        )
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
open(out_path, "r") do f
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
                println("Cycle count: $key")
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

# function to ad N₂ and N₃N₁ vertices and connect edges
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


codon_set =
    LongDNA{4}["CAA", "TTG", "CAC", "GTG", "CAG", "CTG", "CTC", "GAG", "GAA", "TTC", "GAC", "GTC", "GCC"]
data = CodonGraphData(codon_set)
construct_graph_data!(data; show_debug = false)
show_graph(data; show_debug = false)
