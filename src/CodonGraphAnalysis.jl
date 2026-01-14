# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------

"""
    is_c3(data::CodonGraphData; show_debug::Bool = false) -> Bool

Return `true` if the codon set is C3 (original graph and both shifted graphs by 1 and 2 frames are circular).

# Arguments

  - `data::CodonGraphData`: CodonGraphData object containing the codon set.

# Keyword Arguments

  - `show_debug::Bool`: Whether to show debug information (default: false).

# Returns

  - `Bool`: `true` if the codon set is C3, otherwise `false`.

# Throws

  - `ArgumentError`: If the graph has no vertices.
  - `ArgumentError`: If the graph has no edges.
  - `ArgumentError`: If the codon set is empty.
  - `ArgumentError`: If any codon in the codon set is not of length 3.

# Example

```julia
codon_set = LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"])
data = CodonGraphData(codon_set)
construct_graph_data!(data; show_debug = true)    # show original graph
is_c3(data; show_debug = true)
```
"""
function is_c3(data::CodonGraphData; show_debug::Bool = false)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if C3."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if C3."))

    # check if original graph is circular
    if is_circular(data.graph; show_debug = show_debug)
        show_debug && @debug "G(X) is circular -> checking shifted graphs..."

        # check if shifted graphs are circular
        for shift in (1, 2)
            shifted = left_shift_codon_set(data.codon_set, shift)
            shifted_data = CodonGraphData(shifted)
            construct_graph_data!(shifted_data; show_debug = show_debug)
            shift_lower = shift == 1 ? "₁" : "₂"
            if !is_circular(shifted_data.graph; show_debug = show_debug)
                show_debug && @debug "α$shift_lower(X) is not circular -> codon set is not C3"
                return false
            else
                show_debug && @debug "α$shift_lower(X) is circular"
            end
        end
    else
        show_debug && @debug "G(X) is not circular -> codon set is not C3"
        return false
    end
    return true
end


"""
    is_circular(graph::SimpleDiGraph; show_debug::Bool = false) -> Bool

Return `true` if the graph is acyclic (circular).

# Arguments

  - `graph::SimpleDiGraph`: Graph to analyze.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if the graph is acyclic, otherwise `false`.

# Throws

  - `ArgumentError`: If the graph has no vertices.
  - `ArgumentError`: If the graph has no edges.

# Example

```julia
g = SimpleDiGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
is_circular(g)
```
"""
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


"""
    is_graphs_identical(data_1::CodonGraphData, data_2::CodonGraphData; show_debug::Bool = false) -> Bool

Return `true` if both graphs have identical structure and labels.

# Arguments

  - `data_1::CodonGraphData`: First data to compare.
  - `data_2::CodonGraphData`: Second data to compare.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if both graphs match, otherwise `false`.

# Throws

  - None.

# Example

```julia
data_1 = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
data_2 = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data_1)
construct_graph_data!(data_2)
is_graphs_identical(data_1, data_2)
```
"""
function is_codon_graphs_identical(data_1::CodonGraphData, data_2::CodonGraphData; show_debug::Bool = false)
    show_debug && @debug """Comparing graphs...
        Graph 1: nv=$(nv(data_1.graph)), ne=$(ne(data_1.graph))
        Graph 2: nv=$(nv(data_2.graph)), ne=$(ne(data_2.graph))"""
    # check if same amount of vertices and edges
    if nv(data_1.graph) != nv(data_2.graph)
        show_debug && @debug "Not the same amount of vertices"
        return false
    end
    if ne(data_1.graph) != ne(data_2.graph)
        show_debug && @debug "Not the same amount of edges"
        return false
    end

    # check if same vertice labels
    show_debug && @debug "Comparing vertices..."
    for index in 1:nv(data_1.graph)
        show_debug && @debug """In Graph 1: vertice $(index): $(data_1.all_vertex_labels[index])
        In Graph 2: vertice $(index): $(data_2.all_vertex_labels[index])"""
        if !_has_vertex_label(data_2.vertex_index, data_1.all_vertex_labels[index], show_debug = show_debug)
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
        if _has_edge_label(data_2, src_label, dst_label; show_debug = show_debug)
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


"""
    is_comma_free(graph::SimpleDiGraph; show_debug::Bool = false) -> Bool

Return `true` if no path longer than 2 exists in the graph.

# Arguments

  - `graph::SimpleDiGraph`: Graph to analyze.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if the graph is comma-free, otherwise `false`.

# Throws

  - `ArgumentError`: If the graph has no vertices.
  - `ArgumentError`: If the graph has no edges.

# Example

```julia
g = SimpleDiGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
is_comma_free(g)
```
"""
function is_comma_free(graph::SimpleDiGraph; show_debug::Bool = false)
    # do not allow graph with no vertices
    nv(graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if comma-free."))
    # do not allow graphs with no edges
    ne(graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if comma-free."))

    max_depth = 2
    # perform depth-limited DFS for each vertice
    for vertex in vertices(graph)
        if _dfs_depth_limited(graph, vertex, 0, max_depth; show_debug = show_debug)
            show_debug && @debug "Path longer than $max_depth found in graph -> codon set is not comma-free"
            return false
        end
    end

    show_debug && @debug "No paths longer than $max_depth found in graph -> codon set is comma-free"
    return true
end


"""
    is_self_complementary(data::CodonGraphData; show_debug::Bool = false) -> Bool

Return `true` if the graph matches its complemented, reversed graph.

# Arguments

  - `data::CodonGraphData`: Codon graph data to analyze.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if the graph is self-complementary, otherwise `false`.

# Throws

  - `ArgumentError`: If the graph has no vertices.
  - `ArgumentError`: If the graph has no edges.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)    # do not allow graphs with no vertices
is_self_complementary(data)
```
"""
function is_self_complementary(data::CodonGraphData; show_debug::Bool = false)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 &&
        throw(ArgumentError("Graph has no vertices! Cannot determine if self-complementary."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if self-complementary."))


    # create complemented, reversed codon set
    codon_set_complemented_reversed =
        get_complemented_reversed_codon_set(data.codon_set; show_debug = show_debug)

    # create CodonGraphData from complemented, reversed codon set 
    data_complemented_reversed =
        CodonGraphData(codon_set_complemented_reversed; plot_title = "Complemented, reversed graph")
    construct_graph_data!(data_complemented_reversed; show_debug = show_debug)

    # compare original graph with complemented, reversed graph
    if is_codon_graphs_identical(data, data_complemented_reversed; show_debug = show_debug)
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

# function to check if a codon graph is strong C3 (i.e., C3 and expanded graph only has cycles of length 2)
function is_strong_c3(data::CodonGraphData; show_debug::Bool = false)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if strong C3."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if strong C3."))

    # check if original graph is C3
    if is_c3(data; show_debug = show_debug)
        max_depth = 2
        show_debug && @debug "Codon set is C3 -> checking expanded graph for cycles longer than $max_depth..."

        # create expanded graph copy
        data_expanded = CodonGraphData(data.codon_set; plot_title = "Expanded graph for strong C3 check")
        construct_graph_data!(data_expanded; show_debug = show_debug)
        _expand_graph(data_expanded; show_debug = show_debug)

        # check if expanded graph has cycles longer than 2
        for vertex in vertices(data_expanded.graph)
            if _dfs_depth_limited(data_expanded.graph, vertex, 0, max_depth; show_debug = show_debug)
                show_debug &&
                    @debug "Path longer than $max_depth found in expanded graph -> codon set is not strong C3"
                return false
            end
        end

        show_debug &&
            @debug "No paths longer than $max_depth found in expanded graph -> codon set is strong C3"
        return true
    else
        show_debug && @debug "Codon set is not C3 -> codon set cannot be strong C3"
        return false
    end
end


# ---------------------------------------------- HELPERS ----------------------------------------------
# add N₂ and N₃N₁ for each codon N₁N₂N₃ as vertices and connect them accordingly (if possible)
function _add_n2_n3n1_by_codon(data::CodonGraphData, codon::LongDNA{4}; show_debug::Bool = false)
    n2 = string(codon[2])
    n3n1 = string(codon[3], codon[1])
    add_vertex_by_label!(data, n2, show_debug = false)
    add_vertex_by_label!(data, n3n1, show_debug = false)
    add_edge_by_label!(data, n2, n3n1, show_debug = false)
    add_edge_by_label!(data, n3n1, n2, show_debug = false)
    return true
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


# recursive depth-limited DFS to find paths longer than max_depth
function _dfs_depth_limited(
    graph::SimpleDiGraph,
    vertex::Int,
    depth::Int,
    max_depth::Int;
    show_debug::Bool = false,
)
    if depth >= max_depth
        return true
    end

    for neighbor in outneighbors(graph, vertex)
        if _dfs_depth_limited(graph, neighbor, depth + 1, max_depth; show_debug = show_debug)
            return true
        end
    end

    return false
end


# expand graph by adding N₂ and N₃N₁ for each codon N₁N₂N₃ as vertices and connect them accordingly (if possible)
function _expand_graph(data::CodonGraphData; show_debug::Bool = false)
    for codon in data.codon_set
        _add_n2_n3n1_by_codon(data, codon)
    end

    return true
end