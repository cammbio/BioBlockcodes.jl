function is_c3(data::CodonGraphData)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if C3."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if C3."))

    # check if original graph is circular
    if is_circular(data)
        # check if shifted graphs are circular
        for shift in (1, 2)
            shifted_codon_set = left_shift_codon_set(data.codon_set, shift)
            shifted_data = CodonGraphData(shifted_codon_set)
            if !is_circular(shifted_data)
                return false
            end
        end
    else
        return false
    end
    return true
end


function is_circular(data::CodonGraphData)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if circular."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if circular."))

    # state_array to keep track of all vertices (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(data.graph))

    # perform DFS for each vertice
    for vertex in vertices(data.graph)
        if state_array[vertex] == 0
            if _dfs_cycle_detection(data.graph, vertex, state_array)
                return false
            end
        end
    end

    return true
end


function is_comma_free(data::CodonGraphData)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if comma-free."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if comma-free."))

    max_depth = 2
    # perform depth-limited DFS for each vertice
    for vertex in vertices(data.graph)
        if _dfs_depth_limited(data.graph, vertex, 0, max_depth)
            return false
        end
    end

    return true
end


function is_self_complementary(data::CodonGraphData)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 &&
        throw(ArgumentError("Graph has no vertices! Cannot determine if self-complementary."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if self-complementary."))


    # create complemented, reversed codon set
    codon_set_complemented_reversed = get_comp_rev_codon_set(data.codon_set)

    # create CodonGraphData from complemented, reversed codon set 
    data_complemented_reversed = CodonGraphData(codon_set_complemented_reversed)

    # compare original graph with complemented, reversed graph
    if _is_codon_graphs_equal(data, data_complemented_reversed)
        return true
    else
        return false
    end
end


# function to check if a codon graph is strong C3 (i.e., C3 and expanded graph only has cycles of length 2)
function is_strong_c3(data::CodonGraphData)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if strong C3."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if strong C3."))

    # check if original graph is C3
    if is_c3(data)
        max_length = 2

        # create graph copy and expand it
        data_expanded = CodonGraphData(data.codon_set)
        _expand_graph(data_expanded)

        # check if expanded graph has cycles longer than 2
        if _has_cycle_longer_than(data_expanded.graph, max_length)
            return false
        end

        return true
    else
        return false
    end
end


# ---------------------------------------------- HELPERS ----------------------------------------------
# adds edge to graph by vertex labels
function _add_edge_by_label!(data::CodonGraphData, src_label::String, dst_label::String)
    # do not allow empty labels or not existing labels in vert_idxs
    isempty(src_label) && throw(ArgumentError("from_label cannot be empty!"))
    isempty(dst_label) && throw(ArgumentError("to_label cannot be empty!"))
    !haskey(data.vert_idxs, src_label) &&
        throw(ArgumentError("Vertex with label $src_label does not exist in graph!"))
    !haskey(data.vert_idxs, dst_label) &&
        throw(ArgumentError("Vertex with label $dst_label does not exist in graph!"))

    if _has_edge_label(data, src_label, dst_label) # edge already exists
        return false
    else # edge does not already exist
        # get vertex indices from labels
        from_index = data.vert_idxs[src_label]
        to_index = data.vert_idxs[dst_label]

        # add edge to graph and update affected data fields
        push!(data.edge_labels, (src_label, dst_label))
        add_edge!(data.graph, from_index, to_index)
        return true
    end
end


# adds vertices and edges to graph after extracting needed labels from codon set
function _add_vert_by_label!(data::CodonGraphData, label::String)
    # do not allow empty labels
    isempty(label) && throw(ArgumentError("label cannot be empty!"))
    # only allow labels of length 1 or 2
    !(length(label) in (1, 2)) &&
        throw(ArgumentError("label must be of length 1 or 2, got length $(length(label))!"))
    if _has_vert_label(data.vert_idxs, label)
        return false
    else
        add_vertex!(data.graph)
        push!(data.vert_labels, label)
        data.vert_idxs[label] = nv(data.graph)
        return true
    end
end


# recursive DFS to detect cycles
function _dfs_cycle_detection(graph::SimpleDiGraph, vertex::Int, state_array::Vector{Int})
    # mark the current vertice as visiting
    state_array[vertex] = 1
    for neighbor in outneighbors(graph, vertex)
        if state_array[neighbor] == 1
            return true # neighbor already visited, cycle detected
        elseif state_array[neighbor] == 0
            if _dfs_cycle_detection(graph, neighbor, state_array)
                return true # cycle detected in recursion
            end
        end
    end

    # mark the current vertice as visited
    state_array[vertex] = 2
    return false # no cycle detected from this path
end


# recursive DFS to detect cycles longer than min_length
function _dfs_cycles(
    graph::SimpleDiGraph,
    vertex::Int,
    max_length::Int,
    state::Vector{Int},
    stack::Vector{Int},
)
    # mark the current vertice as visiting
    state[vertex] = 1
    # add vertex to current path stack
    push!(stack, vertex)

    for neighbor in outneighbors(graph, vertex)
        if state[neighbor] == 0 # if neighbor not visited, visit it
            _dfs_cycles(graph, neighbor, max_length, state, stack) && return true
        elseif state[neighbor] == 1 # if neighbor is visiting, cycle detected
            index = findlast(==(neighbor), stack)
            cycle_length = length(stack) - index + 1
            if cycle_length > max_length
                return true
            end
        end
    end

    # remove vertex from current path stack
    pop!(stack)
    # mark the current vertice as visited
    state[vertex] = 2
    return false
end


# recursive depth-limited DFS to find paths longer than max_depth
function _dfs_depth_limited(graph::SimpleDiGraph, vertex::Int, depth::Int, max_depth::Int)

    if depth > max_depth
        return true
    end

    # mark the current vertice as visited
    for neighbor in outneighbors(graph, vertex)
        if _dfs_depth_limited(graph, neighbor, depth + 1, max_depth)
            return true
        end
    end

    return false
end


# expand graph by adding N₂ and N₃N₁ for each codon N₁N₂N₃ as vertices and connect them accordingly (if possible)
function _expand_graph(data::CodonGraphData)
    for codon in data.codon_set
        n2 = string(codon[2])
        n3n1 = string(codon[3], codon[1])
        _add_vert_by_label!(data, n2)
        _add_vert_by_label!(data, n3n1)
        _add_edge_by_label!(data, n2, n3n1)
        _add_edge_by_label!(data, n3n1, n2)
        return true
    end

    return true
end


# check if graph has cycle longer than max_length
function _has_cycle_longer_than(graph::SimpleDiGraph, max_length::Int)
    vertex_count = nv(graph)
    on_path = falses(vertex_count)
    stack = Int[]

    function dfs(v::Int, depth::Int)
        push!(stack, v)
        on_path[v] = true

        @inbounds for w in outneighbors(graph, v)
            if on_path[w]
                # w liegt bereits auf dem aktuellen Pfad -> Zyklus gefunden
                idx = findlast(==(w), stack)
                cyc_len = depth - idx + 1
                if cyc_len > max_length
                    pop!(stack)
                    on_path[v] = false
                    return true
                end
                # kürzere Zyklen ignorieren, weiter suchen
            elseif dfs(w, depth + 1)
                pop!(stack)
                on_path[v] = false
                return true
            end
        end

        pop!(stack)
        on_path[v] = false
        return false
    end

    @inbounds for root in 1:vertex_count
        dfs(root, 1) && return true
    end

    return false
end


# check if two codon graphs are equal by comparing vertices and edges
function _is_codon_graphs_equal(data_1::CodonGraphData, data_2::CodonGraphData)
    # check if same amount of vertices and edges
    if nv(data_1.graph) != nv(data_2.graph)
        return false
    end
    if ne(data_1.graph) != ne(data_2.graph)
        return false
    end

    # check if same vertice labels
    for index in 1:nv(data_1.graph)
        if !_has_vert_label(data_2.vert_idxs, data_1.vert_labels[index])
            return false
        end
    end

    # check if same edges
    for edge in edges(data_1.graph)
        src_label = data_1.vert_labels[src(edge)]
        dst_label = data_1.vert_labels[dst(edge)]
        if _has_edge_label(data_2, src_label, dst_label)
        else
            return false
        end
    end
    return true
end


# checks if an edge with given labels exists in the graph
function _has_edge_label(data::CodonGraphData, src_label::String, dst_label::String)
    # check if labels exist
    haskey(data.vert_idxs, src_label) || return false
    haskey(data.vert_idxs, dst_label) || return false

    src_index = data.vert_idxs[src_label]
    dst_index = data.vert_idxs[dst_label]

    return has_edge(data.graph, src_index, dst_index)
end


# checks if a vertex with given label exists in the graph
function _has_vert_label(vert_idxs::Dict{String, Int}, label::String)
    return haskey(vert_idxs, label)
end
