# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
# function to get count of all cycles in a graph
function get_cycle_count(graph::SimpleDiGraph; debug::Bool = false)
    cycles = simplecycles(graph)
    return length(cycles)
end


# function to get count of cycles with length equal to cycle_length
function get_cycle_count_by_length(data::CodonGraphData, cycle_length::Int; debug::Bool = false)
    cycles = simplecycles(data.graph)
    cycles_filtered = filter(cycle -> length(cycle) == cycle_length, cycles)
    return length(cycles_filtered)
end


# function to get all cycles in a graph
function get_cycles_all(data::CodonGraphData; debug::Bool = false)
    # get all cycles using simplecycles
    cycles = simplecycles(data.graph)

    # create new array to store cycles with vertice labels instead of indices
    cycles_labeled = Vector{Vector{String}}()
    # iterate all cycles and convert vertice indices to labels
    for cycle in cycles
        cycle_labeled = [data.all_vertex_labels[i] for i in cycle]
        push!(cycles_labeled, cycle_labeled)
    end
    return cycles_labeled
end


# function to get all cycles with length equal to cycle_length
function get_cycles_by_length(data::CodonGraphData, cycle_length::Int; debug::Bool = false)
    cycles = simplecycles(data.graph)
    cycles_filtered = filter(cycle -> length(cycle) == cycle_length, cycles)
    # iterate all filtered cycles and convert vertice indices to labels
    cycles_labeled = Vector{Vector{String}}()
    for cycle in cycles_filtered
        cycle_labeled = [data.all_vertex_labels[i] for i in cycle]
        push!(cycles_labeled, cycle_labeled)
    end
    return cycles_labeled
end


# function to get all cycles that start from a given vertex label in a graph
function get_cycles_by_vertex_label(data::CodonGraphData, vertex_label::String; debug::Bool = false)
    !haskey(data.vertex_index, vertex_label) &&
        throw(ArgumentError("Vertex label $vertex_label not found in graph data."))

    vertex_index = data.vertex_index[vertex_label]
    cycles = simplecycles(data.graph)
    cycles_filtered = filter(cycle -> cycle[1] == vertex_index, cycles)
    return cycles_filtered
end


# function to get cycles which are in graph1 but not in graph2
function get_cycles_difference(data_1::CodonGraphData, data_2::CodonGraphData; debug::Bool = false)
    # get cycles from both graphs
    cycles_1 = simplecycles(data_1.graph)
    cycles_2 = simplecycles(data_2.graph)

    # convert cycles to sets of tuples for easy comparison
    cycles_set_1 = Set(map(Tuple, cycles_1))
    cycles_set_2 = Set(map(Tuple, cycles_2))

    # get difference of cycles
    difference_keys = setdiff(cycles_set_1, cycles_set_2)
    return difference_keys
end


# function to display duplicate cycles
function get_duplicate_cycles(data::CodonGraphData; debug::Bool = false)
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


# function to get the maximal cycle length in the graph
function get_max_cycle_length(graph::SimpleDiGraph; debug::Bool = false)
    cycles = simplecycles(graph)
    isempty(cycles) && return 0
    return maximum(length.(cycles))
end
