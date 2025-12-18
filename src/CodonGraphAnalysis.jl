# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
is_circular(data::CodonGraphData) -> Bool

Returns true if the graph is circular aka. acyclic
"""
# check if a set of codons is circular by checking if the graph is acyclic
function is_circular(graph::SimpleDiGraph; show_debug::Bool = false)
    # do not allow graph with no vertices
    nv(graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if circular."))
    # do not allow graphs with no edges
    ne(graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if circular."))
    # state_array to keep track of all vertices (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(graph))

    # perform DFS for each vertice
    for vertex in vertices(graph)
        if state_array[vertex] == 0
            if _dfs_cycle_detection(graph, vertex, state_array; show_debug = show_debug)
                show_debug && @debug "Cycle detected in graph -> Graph is not acyclic aka. not circular"
                return false # no cycles detected, Graph is acyclic aka. circular
            end
        end
    end

    show_debug && @debug "No cycles detected in graph -> Graph is acyclic aka. circular"
    return true # no cycles detected, Graph is not acyclic aka. not circular
end


# check if a set of codons is C3 by checking if the two shifted graphs α₁(X) and α₂ are circular
function is_c3(data::CodonGraphData; show_plot::Bool = false, show_debug::Bool = false)
    # show original graph
    if show_plot
        show_graph(data; show_debug = show_debug)
    end
    # create shifted graph
    shifted_data_by_1 =
        create_shifted_graph(data.codon_set, 1; show_debug = show_debug, show_plot = show_plot)
    shifted_data_by_2 =
        create_shifted_graph(data.codon_set, 2; show_debug = show_debug, show_plot = show_plot)

    # check if original graph and both shifted graphs are circular
    if is_circular(data.graph; show_debug = show_debug) &&
       is_circular(shifted_data_by_1.graph; show_debug = show_debug) &&
       is_circular(shifted_data_by_2.graph; show_debug = show_debug)
        show_debug && @debug "G(X), α₁(X) and α₂ are circular -> codon set is C3"
        return true
    else
        show_debug && @debug "G(X), α₁(X) or α₂(X) is not circular -> codon set is not C3"
        return false
    end
end


# check if a set of codons is comma-free by checking if a path longer than 2 exists
function is_comma_free(graph::SimpleDiGraph; show_debug::Bool = false)
    # do not allow graph with no vertices
    nv(graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if comma-free."))
    # do not allow graphs with no edges
    ne(graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if comma-free."))

    # perform depth-limited DFS for each vertice
    for vertex in vertices(graph)
        if _dfs_depth_limited(graph, vertex, 0; show_debug = show_debug)
            show_debug && @debug "Path longer than 2 found in graph -> codon set is not comma-free"
            return false # path longer than 2 found
        end
    end

    show_debug && @debug "No paths longer than 2 found in graph -> codon set is comma-free"
    return true # no paths longer than 2 found
end




# check if a set of codons is self-complementary by checking if the graph G(X) is identical to its complemented, reversed graph
function is_self_complementary(data::CodonGraphData; show_debug::Bool = false, show_plot::Bool = false)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 &&
        throw(ArgumentError("Graph has no vertices! Cannot determine if self-complementary."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if self-complementary."))

    # create complemented, reversed codon set
    codon_set_complemented_reversed =
        get_complemented_reversed_codon_set(data.codon_set; show_debug = show_debug)

    # create complemented, reversed graph
    data_complemented_reversed =
        CodonGraphData(codon_set_complemented_reversed; plot_title = "Complemented, reversed graph")
    construct_graph!(data_complemented_reversed; show_debug = show_debug, show_plot = show_plot)

    # compare original graph with complemented, reversed graph
    if is_graphs_identical(data, data_complemented_reversed; show_debug = show_debug)
        show_debug && @debug """Original codon set:
        $(data.codon_set)
        Complemented, reversed codon set:
        $(data_complemented_reversed.codon_set)
        Graphs are identical -> original codon set is self-complementary"""
        return true
    else
        show_debug && @debug """Original codon set:
        $(data.codon_set)
        Complemented, reversed codon set:
        $(data_complemented_reversed.codon_set)
        Graphs are not identical -> original codon set is not self-complementary"""
        return false
    end
end


# compare two graphs for identical structure
function is_graphs_identical(data_1::CodonGraphData, data_2::CodonGraphData; show_debug::Bool = false)
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
        if !has_vertex_label(data_2.vertex_index, data_1.all_vertex_labels[index], show_debug = show_debug)
            show_debug && @debug "vertice label $(data_1.all_vertex_labels[index]) NOT found in Graph 2"
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


# check if a vertex label exists in the graph
function has_vertex_label(vertex_index::Dict{String, Int}, label::String; show_debug::Bool = false)
    show_debug && @debug "Checking if vertice label $label exists in graph..."
    return haskey(vertex_index, label)
end


# check if edge exists in graph between two vertice labels
function has_edge_label(data::CodonGraphData, src_label::String, dst_label::String; show_debug::Bool = false)
    # check if labels exist
    haskey(data.vertex_index, src_label) || return false
    haskey(data.vertex_index, dst_label) || return false

    from_index = data.vertex_index[src_label]
    to_index = data.vertex_index[dst_label]

    return has_edge(data.graph, from_index, to_index)
end


# recursive DFS to detect cycles
function _dfs_cycle_detection(
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
            if _dfs_cycle_detection(graph, neighbor, state_array; show_debug = show_debug)
                return true # cycle detected in recursion
            end
        end
    end

    # mark the current vertice as visited
    state_array[vertex] = 2
    return false # no cycle detected from this path
end


# recursive depth-limited DFS to find paths longer than 2
function _dfs_depth_limited(graph::SimpleDiGraph, vertex::Int, depth::Int; show_debug::Bool = false)
    if depth >= 3
        return true # Pfad mit Länge >= 3 gefunden
    end

    for neighbor in outneighbors(graph, vertex)
        if _dfs_depth_limited(graph, neighbor, depth + 1; show_debug = show_debug)
            return true
        end
    end
    return false
end