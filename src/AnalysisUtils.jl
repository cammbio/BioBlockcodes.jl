# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
"""
is_circular(data::CodonGraphData) -> Bool

Returns true if the codon graph represented by `data` is circular aka. is acyclic (does not contain any cycles)
"""
# check if a set of codons is circular by checking if the graph is acyclic
function is_circular(graph::SimpleDiGraph; show_debug::Bool = false)
    # state_array to keep track of all vertices (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(graph))

    # perform DFS for each vertice
    for vertex in vertices(graph)
        if state_array[vertex] == 0
            if dfs_cycle_detection(graph, vertex, state_array; show_debug = show_debug)
                println("Cycle detected in graph -> Graph is not acyclic aka. not circular")
                return false # no cycles detected, Graph is acyclic aka. circular
            end
        end
    end

    println("No cycles detected in graph -> Graph is acyclic aka. circular")
    return true # no cycles detected, Graph is not acyclic aka. not circular
end

# recursive DFS to detect cycles
function dfs_cycle_detection(
    graph::SimpleDiGraph,
    vertex::Int,
    state_array::Vector{Int};
    show_debug::Bool = false,
)
    # mark the current vertice as visiting
    state_array[vertex] = 1
    for neighbor in outneighbors(graph, vertex)
        if state_array[neighbor] == 1
            return true # neighbor already visited, cycle detected
        elseif state_array[neighbor] == 0
            if dfs_cycle_detection(graph, neighbor, state_array; show_debug = show_debug)
                return true # cycle detected in recursion
            end
        end
    end

    # mark the current vertice as visited
    state_array[vertex] = 2
    return false # no cycle detected from this path
end


# check if a set of codons is comma-free by checking if a path longer than 2 exists
function is_comma_free(graph::SimpleDiGraph; show_debug::Bool = false)
    for vertex in vertices(graph)
        if dfs_depth_limited(graph, vertex, 0; show_debug = show_debug)
            println("Path longer than 2 found in graph -> codon set is not comma-free")
            return false # path longer than 2 found
        end
    end

    println("No paths longer than 2 found in graph -> codon set is comma-free")
    return true # no paths longer than 2 found
end


# recursive depth-limited DFS to find paths longer than 2
function dfs_depth_limited(graph::SimpleDiGraph, vertex::Int, depth::Int; show_debug::Bool = false)
    if depth >= 3
        return true # Pfad mit Länge >= 3 gefunden
    end

    for neighbor in outneighbors(graph, vertex)
        if dfs_depth_limited(graph, neighbor, depth + 1; show_debug = show_debug)
            return true
        end
    end

    return false
end


# check if a set of codons is self-complementary by checking if the graph G(X) is equal to its inverted Graph
function is_self_complementary(
    data::CodonGraphData;
    show_plot::Bool = true,
    show_debug::Bool = false,
)
    # create complement reversed codons
    codon_set_complemented_reversed =
        get_complemented_reversed_codon_set(data; show_debug = show_debug)
    # create new graph with complement_reversed_codon_set to compare with original graph
    data_complemented_reversed =
        CodonGraphData(codon_set_complemented_reversed; plot_title = "Complemented, reversed graph")
    construct_graph!(data_complemented_reversed; show_plot = show_plot, show_debug = show_debug)
    if is_graphs_identical(data, data_complemented_reversed; show_debug = show_debug)
        println("""Original codon set:
        $(data.codon_set)
        Complemented, reversed codon set:
        $(data_complemented_reversed.codon_set)
        Graphs are identical -> original codon set is self-complementary""")
        return true
    else
        println("""Original codon set:
        $(data.codon_set)
        Complemented, reversed codon set:
        $(data_complemented_reversed.codon_set)
        Graphs are not identical -> original codon set is not self-complementary""")
        return false
    end
end


# compare two graphs for identical structure
function is_graphs_identical(
    data_1::CodonGraphData,
    data_2::CodonGraphData;
    show_debug::Bool = false,
)
    # check if same amount of vertices and edges
    show_debug && @debug """Comparing graphs...
    Graph 1: nv=$(nv(data_1.graph)), ne=$(ne(data_1.graph))
    Graph 2: nv=$(nv(data_2.graph)), ne=$(ne(data_2.graph))"""
    if nv(data_1.graph) != nv(data_2.graph)
        show_debug && @debug "Not the same amount of vertices"
        return false
    end
    if ne(data_1.graph) != ne(data_2.graph)
        show_debug && @debug "Not the same amount of edges"
        return false
    end

    # check if same vertice labels
    show_debug && @debug "Comparing vertice labels..."
    for index in 1:nv(data_1.graph)
        show_debug && @debug """In Graph 1: vertice $(index): $(data_1.all_vertex_labels[index])
        In Graph 2: vertice $(index): $(data_2.all_vertex_labels[index])"""
        if !has_vertex_label(
            data_2.vertex_index,
            data_1.all_vertex_labels[index],
            show_debug = show_debug,
        )
            show_debug &&
                @debug "vertice label $(data_1.all_vertex_labels[index]) NOT found in Graph 2"
            return false
        end
    end

    # check if same edges
    show_debug && @debug "Comparing edges..."
    for edge in edges(data_1.graph)
        src_label = data_1.all_vertex_labels[src(edge)]
        dst_label = data_1.all_vertex_labels[dst(edge)]
        show_debug &&
            @debug "Edge: $(data_1.all_vertex_labels[src(edge)]) -> $(data_1.all_vertex_labels[dst(edge)])"
        if has_edge_label(data_2, src_label, dst_label; show_debug = show_debug)
            show_debug &&
                @debug "Edge: $(data_1.all_vertex_labels[src(edge)]) -> $(data_1.all_vertex_labels[dst(edge)]) also in Graph 2"
        else
            show_debug &&
                @debug "Edge: $(data_1.all_vertex_labels[src(edge)]) -> $(data_1.all_vertex_labels[dst(edge)]) NOT in Graph 2"
            return false # edge not found
        end
    end

    return true # graphs are identical
end


# check if edge exists in graph between two vertice labels
function has_edge_label(
    data::CodonGraphData,
    src_label::String,
    dst_label::String;
    show_debug::Bool = false,
)
    # check if labels exist
    haskey(data.vertex_index, src_label) || return false
    haskey(data.vertex_index, dst_label) || return false

    from_index = data.vertex_index[src_label]
    to_index = data.vertex_index[dst_label]

    return has_edge(data.graph, from_index, to_index)
end


# check if a vertice labels exists in the graph
function has_vertex_label(vertex_index::Dict{String, Int}, label::String; show_debug::Bool = false)
    return haskey(vertex_index, label)
end


# check if a set of codons is C3 by checking if the two shifted graphs α₁(X) and α₂ are circular
function is_c3(data::CodonGraphData; show_plot::Bool = false, show_debug::Bool = false)
    # show original graph
    if show_plot
        show_graph(data; show_debug = show_debug)
    end
    # create shifted graph
    shifted_data_by_1 =
        create_shifted_graph(data.codon_set, 1; show_plot = show_plot, show_debug = show_debug)
    shifted_data_by_2 =
        create_shifted_graph(data.codon_set, 2; show_plot = show_plot, show_debug = show_debug)

    # check if original graph and both shifted graphs are circular
    if is_circular(data.graph; show_debug = show_debug) &&
       is_circular(shifted_data_by_1.graph; show_debug = show_debug) &&
       is_circular(shifted_data_by_2.graph; show_debug = show_debug)
        println("G(X), α₁(X) and α₂ are circular -> codon set is C3")
        return true
    else
        println("G(X), α₁(X) or α₂(X) is not circular -> codon set is not C3")
        return false
    end
end


# create shifted graph αₖ(X) from original graph by shifting codons by k positions
function create_shifted_graph(
    codon_set::Vector{LongDNA{4}},
    shift_by::Int;
    show_plot::Bool = false,
    show_debug::Bool = false,
)
    # create shifted codon set
    shifted_codon_set = left_shift_codon_set(codon_set, shift_by; show_debug = show_debug)

    # create new CodonGraphData for shifted graph
    shifted_data = CodonGraphData(shifted_codon_set; plot_title = "Shifted graph by $shift_by")
    construct_graph!(shifted_data; show_plot = show_plot, show_debug = show_debug)
    return shifted_data
end


# shift a codon set by k positions to the left
function left_shift_codon_set(
    codon_set::Vector{LongDNA{4}},
    shift_by::Int;
    show_debug::Bool = false,
)
    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon_set[1]))
    # shift every codon from codon_set
    shifted_codon_set = Vector{LongDNA{4}}()
    for codon in codon_set
        # cut of first shift_by characters and append them to the end
        shifted_codon = left_shift_codon(codon, shift_by; show_debug = show_debug)
        push!(shifted_codon_set, shifted_codon)
    end
    show_debug && @debug """Original codon set: $codon_set
    -> shifted codon set by $shift_by: $shifted_codon_set"""
    return shifted_codon_set
end


# shift a codon by k positions to the left
function left_shift_codon(codon::LongDNA{4}, shift_by::Int; show_debug::Bool = false)
    # limit shift_by to length of codon
    shift_by = mod(shift_by, length(codon))
    # cut of first shift_by characters and append them to the end
    shifted_codon = codon[(shift_by + 1):end] * codon[1:shift_by]
    show_debug && @debug """Original codon: $codon
    -> shifted codon by $shift_by: $shifted_codon"""
    return shifted_codon
end


# function to add a new vertice to a graph data structure
function add_vertex_by_label!(data::CodonGraphData, label::String; show_debug::Bool = false)
    if label in data.all_vertex_labels # vertice already exists
        show_debug && @debug "Vertice $label already exists in graph -> not added."
        return false
    else # vertice does not already exist
        # update affected data fields
        add_vertex!(data.graph) # add vertice to graph
        push!(data.added_vertice_labels, label) # add to manually added vertice labels
        data.vertex_index[label] = nv(data.graph) # map label to vertice index
        show_debug && @debug "Added vertice: $label"
        return true
    end
end


# function to add a new edge to a graph data structure
function add_edge_by_label!(
    data::CodonGraphData,
    from_label::String,
    to_label::String;
    show_debug::Bool = false,
)
    if has_edge_label(data, from_label, to_label; show_debug = show_debug) # edge already exists
        show_debug && @debug "Edge $from_label -> $to_label already exists in graph -> not added."
        return false
    else # edge does not already exist
        connect_edge_by_label!(data, from_label, to_label; show_debug = show_debug)
        push!(data.added_edge_labels, (from_label, to_label)) # add to manually added edge labels
        show_debug && @debug "Added edge: $from_label -> $to_label"
        return true
    end
end


# connect one edge label to another
function connect_edge_by_label!(
    data::CodonGraphData,
    from_label::String,
    to_label::String;
    show_debug::Bool = false,
)
    # get needed vertice IDs
    from_index = data.vertex_index[from_label]
    to_index = data.vertex_index[to_label]
    add_edge!(data.graph, from_index, to_index)
end


# show all cycles in the graph
function display_cycles(data::CodonGraphData; show_debug::Bool = false)
    cycles = simplecycles(data.graph)
    show_debug && @debug """Found $(length(cycles)) cycles in graph.
    Cycles: $cycles"""

    # check for duplicates
    keys = map(Tuple, cycles)
    unique_keys = unique(keys)
    show_debug && @debug """Unique cycles: $unique_keys
    Total cycles: $keys"""
    if length(keys) != length(unique_keys)
        show_debug && @debug "Duplicate cycles found!"
        return false
    end

    # iterate all cycles and print them
    for cycle in cycles # iterate all cycles
        # get all vertice labels 
        # join every vertice label in the cycle with " -> " and print it
        println(join((data.all_vertex_labels[i] for i in (cycle..., first(cycle))), " -> "))
        println("Cycle length: $(length(cycle))")
    end
    println("Amount of cycles found: $(length(cycles))")
end