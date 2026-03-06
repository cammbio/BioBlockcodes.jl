# DEPRECATED
# function to get count of all cycles in a graph
function get_cycle_count(graph::SimpleDiGraph)
    cycles = simplecycles(graph)
    return length(cycles)
end


# function to get count of cycles with length equal to cycle_length
function get_cycle_count_by_length(cgd::CodonGraphData, cycle_length::Int)
    cycles = simplecycles(cgd.graph)
    cycles_filtered = filter(cycle -> length(cycle) == cycle_length, cycles)
    return length(cycles_filtered)
end


# function to get all cycles in a graph
function get_cycles_all(cgd::CodonGraphData)
    # get all cycles using simplecycles
    cycles = simplecycles(cgd.graph)

    # create new array to store cycles with vertice labels instead of indices
    cycles_labeled = Vector{Vector{String}}()
    # iterate all cycles and convert vertice indices to labels
    for cycle in cycles
        cycle_labeled = [cgd.vert_labels[i] for i in cycle]
        push!(cycles_labeled, cycle_labeled)
    end
    return cycles_labeled
end


# function to get all cycles with length equal to cycle_length
function get_cycles_by_length(cgd::CodonGraphData, cycle_length::Int)
    cycles = simplecycles(cgd.graph)
    cycles_filtered = filter(cycle -> length(cycle) == cycle_length, cycles)
    # iterate all filtered cycles and convert vertice indices to labels
    cycles_labeled = Vector{Vector{String}}()
    for cycle in cycles_filtered
        cycle_labeled = [cgd.vert_labels[i] for i in cycle]
        push!(cycles_labeled, cycle_labeled)
    end
    return cycles_labeled
end


# function to get all cycles that start from a given vertex label in a graph
function get_cycles_by_vertex_label(cgd::CodonGraphData, vertex_label::String)
    !haskey(cgd.vert_idxs, vertex_label) &&
        throw(ArgumentError("Vertex label $vertex_label not found in graph cgd."))

    vert_idxs = cgd.vert_idxs[vertex_label]
    cycles = simplecycles(cgd.graph)
    cycles_filtered = filter(cycle -> cycle[1] == vert_idxs, cycles)
    return cycles_filtered
end


# function to get cycles which are in graph1 but not in graph2
function get_cycles_difference(cgd_1::CodonGraphData, cgd_2::CodonGraphData)
    # get cycles from both graphs
    cycles_1 = simplecycles(cgd_1.graph)
    cycles_2 = simplecycles(cgd_2.graph)

    # convert cycles to sets of tuples for easy comparison
    cycles_set_1 = Set(map(Tuple, cycles_1))
    cycles_set_2 = Set(map(Tuple, cycles_2))

    # get difference of cycles
    difference_keys = setdiff(cycles_set_1, cycles_set_2)
    return difference_keys
end


# function to display duplicate cycles
function get_duplicate_cycles(cgd::CodonGraphData)
    println("Checking for duplicate cycles in the graph...")
    cycles = simplecycles(cgd.graph)
    cycle_keys = map(Tuple, cycles)
    seen = Set{Tuple{Vararg{Int64}}}()
    duplicate_keys = Tuple{Vararg{Int64}}[]
    for key in cycle_keys
        if key in seen
            push!(duplicate_keys, key)
            println(
                "Duplicate cycle found: ",
                join((cgd.vert_labels[i] for i in (key..., first(key))), " -> "),
            )
        else
            push!(seen, key)
        end
    end
end


# function to get the maximal cycle length in the graph
function get_max_cycle_length(graph::SimpleDiGraph)
    cycles = simplecycles(graph)
    isempty(cycles) && return 0
    return maximum(length.(cycles))
end
