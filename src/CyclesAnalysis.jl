# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# function to display all cycles in the graph
function display_all_cycles(data::CodonGraphData; show_debug::Bool = false)
    cycles = simplecycles(data.graph)

    # iterate all cycles and print them
    for cycle in cycles # iterate all cycles
        # join every vertice label in the cycle with " -> " and print it
        println(join((data.all_vertex_labels[i] for i in (cycle..., first(cycle))), " -> "))
        println("Cycle length: $(length(cycle))")
    end
end


# function to display duplicate cycles
function display_duplicate_cycles(data::CodonGraphData; show_debug::Bool = false)
    println("Checking for duplicate cycles in the graph...")
    cycles = simplecycles(data.graph)
    cycle_keys = map(Tuple, cycles)
    seen = Set{Tuple{Vararg{Int64}}}()
    duplicate_keys = Tuple{Vararg{Int64}}[]
    for key in cycle_keys
        if key in seen
            push!(duplicate_keys, key)
            println(
                "Duplicate cycle found: ",
                join((data.all_vertex_labels[i] for i in (key..., first(key))), " -> "),
            )
        else
            push!(seen, key)
        end
    end
end


# function to get all cycles that start from a given vertex
function get_cycles_from_vertex(data::CodonGraphData, vertex_label::String; show_debug::Bool = false)
    !haskey(data.vertex_index, vertex_label) &&
        throw(ArgumentError("Vertex label $vertex_label not found in graph data."))

    vertex_index = data.vertex_index[vertex_label]
    cycles = simplecycles(data.graph)
    filtered_cycles = filter(cycle -> cycle[1] == vertex_index, cycles)
    println("Cycles starting from vertex $vertex_label:")
    for cycle in filtered_cycles
        println(join((data.all_vertex_labels[i] for i in (cycle..., first(cycle))), " -> "))
        # println("Cycle length: $(length(cycle))")
    end
    return filtered_cycles
end


# function to get count of all cycles
function get_cycle_count(graph::SimpleDiGraph; show_debug::Bool = false)
    cycles = simplecycles(graph)
    return length(cycles)
end


# function to get all cycles with length equal to cycle_length
function get_cycles_of_length(data::CodonGraphData, cycle_length::Int; show_debug::Bool = false)
    cycles = simplecycles(data.graph)
    filtered_cycles = filter(cycle -> length(cycle) == cycle_length, cycles)
    println("Cycles with length $cycle_length:")
    for cycle in filtered_cycles
        println(join((data.all_vertex_labels[i] for i in (cycle..., first(cycle))), " -> "))
    end
end


# function to get count of cycles with length equal to cycle_length
function get_cycle_count_of_length(data::CodonGraphData, cycle_length::Int; show_debug::Bool = false)
    cycles = simplecycles(data.graph)
    filtered_cycles = filter(cycle -> length(cycle) == cycle_length, cycles)
    return length(filtered_cycles)
end


# function to get cycles between two graphs which are only in one of the graphs,
# e.g. in graph1 but not in graph2 and vice versa
function get_cycles_difference(data1::CodonGraphData, data2::CodonGraphData; show_debug::Bool = false)
    cycles1 = simplecycles(data1.graph)
    cycles2 = simplecycles(data2.graph)

    cycle_keys1 = Set(map(Tuple, cycles1))
    cycle_keys2 = Set(map(Tuple, cycles2))

    difference_keys1 = setdiff(cycle_keys1, cycle_keys2)
    difference_keys2 = setdiff(cycle_keys2, cycle_keys1)

    println("Cycles in graph1 which not in graph2:")
    for key in difference_keys1
        println(join((data1.all_vertex_labels[i] for i in (key..., first(key))), " -> "))
    end

    println("Cylces in graph2 which not in graph1:")
    for key in difference_keys2
        println(join((data2.all_vertex_labels[i] for i in (key..., first(key))), " -> "))
    end
end


# function to check cycles of all lengths in the graph
function check_cycle_count_by_length(
    data::CodonGraphData;
    start_index::Int = 1,
    end_index::Int = 100,
    step_size::Int = 1,
    show_debug::Bool = false,
)
    println("Cycles count: $(get_cycle_count(data; show_debug = false))")
    for i in start_index:step_size:end_index
        if get_cycle_count_of_length(data, i; show_debug = show_debug) == 0
            continue
        end
        println("Count of cycles of length $i: $(get_cycle_count_of_length(data, i; show_debug = true))")
    end
end


# function to get the maximal cycle length in the graph
function get_max_cycle_length(graph::SimpleDiGraph; show_debug::Bool = false)
    cycles = simplecycles(graph)
    max_length = maximum(length.(cycles))
    return max_length
end