function is_c3(data::CodonGraphData)
    # validate data object
    _validate_cgd(data)
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
    # validate data object
    _validate_cgd(data)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if circular."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if circular."))

    # state_array to keep track of all vertices (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(data.graph))

    # perform DFS for each vertice
    for vertex in vertices(data.graph)
        if state_array[vertex] == 0
            if _has_cycle_longer_than(data.graph, 0)
                return false
            end
        end
    end

    return true
end


function is_comma_free(data::CodonGraphData)
    # validate data object
    _validate_cgd(data)
    # do not allow graph with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if comma-free."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if comma-free."))

    max_depth = 2
    # perform depth-limited DFS for each vertice
    for vertex in vertices(data.graph)
        if _dfs_depth(data.graph, vertex, 0, max_depth)
            return false
        end
    end

    return true
end


function is_self_complementary(data::CodonGraphData)
    # validate data object
    _validate_cgd(data)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 &&
        throw(ArgumentError("Graph has no vertices! Cannot determine if self-complementary."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if self-complementary."))

    # create complemented, reversed codon set
    comp_rev_codon_set = get_comp_rev_codon_set(data.codon_set)

    # create CodonGraphData from complemented, reversed codon set 
    comp_rev_data = CodonGraphData(comp_rev_codon_set)
    _validate_cgd(comp_rev_data)

    # compare original graph with complemented, reversed graph
    if _is_codon_graphs_equal(data, comp_rev_data)
        return true
    else
        return false
    end
end


# function to check if a codon graph is strong C3 (i.e., C3 and expanded graph only has cycles of length 2)
function is_strong_c3(data::CodonGraphData)
    # validate data object
    _validate_cgd(data)
    # do not allow graphs with no vertices
    nv(data.graph) == 0 && throw(ArgumentError("Graph has no vertices! Cannot determine if strong C3."))
    # do not allow graphs with no edges
    ne(data.graph) == 0 && throw(ArgumentError("Graph has no edges! Cannot determine if strong C3."))

    # check if original graph is C3
    if is_c3(data)
        # create graph copy
        data_expanded = CodonGraphData(data.codon_set)
        _validate_cgd(data_expanded)
        # expand copy graph
        _expand_graph(data_expanded)
        # validate expanded data object
        _validate_cgd(data_expanded)

        # check if expanded graph has cycles longer than 2
        max_length = 2
        if _has_cycle_longer_than(data_expanded.graph, max_length)
            return false
        end

        return true
    else
        return false
    end
end


# adds edge to graph by vertex labels
function _add_edge_by_label!(data::CodonGraphData, src_label::String, dst_label::String)
    # do not allow empty labels
    isempty(src_label) && throw(ArgumentError("2nd argument cannot be empty!"))
    isempty(dst_label) && throw(ArgumentError("3rd argument cannot be empty!"))
    # labels must be of length 1 or 2
    !(length(src_label) in (1, 2)) &&
        throw(ArgumentError("2nd argument must be of length 1 or 2, got length $(length(src_label))!"))
    !(length(dst_label) in (1, 2)) &&
        throw(ArgumentError("3rd argument must be of length 1 or 2, got length $(length(dst_label))!"))
    # check if labels exist in vert_idxs
    !haskey(data.vert_idxs, src_label) &&
        throw(ArgumentError("Vertex with label $src_label does not exist in graph!"))
    !haskey(data.vert_idxs, dst_label) &&
        throw(ArgumentError("Vertex with label $dst_label does not exist in graph!"))


    # edge already exists
    if _has_edge_label(data, src_label, dst_label)
        println("Edge from $src_label to $dst_label already exists in graph! No edge added.")
        return false
    else
        src_idx = data.vert_idxs[src_label]
        dst_idx = data.vert_idxs[dst_label]
        # add edge to graph and update affected data fields
        if add_edge!(data.graph, src_idx, dst_idx)
            push!(data.edge_labels, (src_label, dst_label))
            println("Edge from $src_label to $dst_label added to graph.")
            return true
        end
    end
end


# adds vertices and edges to graph after extracting needed labels from codon set
function _add_vert_by_label!(data::CodonGraphData, label::String)
    # do not allow empty labels
    isempty(label) && throw(ArgumentError("label cannot be empty."))
    # only allow labels of length 1 or 2
    !(length(label) in (1, 2)) &&
        throw(ArgumentError("label must be of length 1 or 2, got length $(length(label))."))
    # only allow labels containing A, C, G or T
    for char in label
        char in ALLOWED_BASES_STR || throw(
            ArgumentError("label contains invalid character \"$char\". Only A, C, G and T are allowed."),
        )
    end

    if _has_vert_label(data.vert_idxs, label)
        return false
    else
        add_vertex!(data.graph)
        push!(data.vert_labels, label)
        data.vert_idxs[label] = nv(data.graph)
        return true
    end
end


# recursive depth-limited DFS to find paths longer than max_depth
function _dfs_depth(graph::SimpleDiGraph, vertex::Int, depth::Int, max_depth::Int)
    # do not allow negative max_depth
    max_depth < 0 && throw(ArgumentError("max_depth must be at least 0."))

    if depth > max_depth
        return true
    end

    # mark the current vertice as visited
    for neighbor in outneighbors(graph, vertex)
        # recursively visit neighbor with increased depth
        if _dfs_depth(graph, neighbor, depth + 1, max_depth)
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
    end
end


# check if graph has cycle longer than max_length
function _has_cycle_longer_than(graph::SimpleDiGraph, max_length::Int)
    max_length < 0 && throw(ArgumentError("max_length must be at least 0"))

    # stack for DFS: (current_vertex, remaining_neighbors, current_depth)
    dfs_stack = Vector{Tuple{Int, Vector{Int}, Int}}()
    vert_count = nv(graph)
    curr_path = Int[]
    in_curr_path = falses(vert_count)

    for root_idx in 1:vert_count
        # reset stack and path tracking for each new root vertex
        empty!(dfs_stack)
        empty!(curr_path)
        fill!(in_curr_path, false)

        # initialize stack with root vertex
        push!(dfs_stack, (root_idx, collect(outneighbors(graph, root_idx)), 1))
        push!(curr_path, root_idx)
        in_curr_path[root_idx] = true
        # DFS loop
        while !isempty(dfs_stack)
            curr_vert, rem_neighs, curr_depth = dfs_stack[end]

            # if no more neighbors to visit, backtrack
            if isempty(rem_neighs)
                # backtrack
                pop!(dfs_stack)
                pop!(curr_path)
                in_curr_path[curr_vert] = false
                continue
            end

            # visit next neighbor
            next_neigh = popfirst!(rem_neighs)
            dfs_stack[end] = (curr_vert, rem_neighs, curr_depth)

            if in_curr_path[next_neigh]
                idx = findlast(==(next_neigh), curr_path)
                cyc_len = curr_depth - idx + 1
                cyc_len > max_length && return true
            else
                push!(dfs_stack, (next_neigh, collect(outneighbors(graph, next_neigh)), curr_depth + 1))
                push!(curr_path, next_neigh)
                in_curr_path[next_neigh] = true
            end
        end
    end

    return false
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


# check if two codon graphs are equal by comparing vertices and edges
function _is_codon_graphs_equal(data_1::CodonGraphData, data_2::CodonGraphData)
    # check if same amount of vertices
    if nv(data_1.graph) != nv(data_2.graph)
        return false
    end
    # check if same amount of edges
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


# validate CodonGraphData consistency (labels, indices, graph size)
function _validate_cgd(data::CodonGraphData)
    _validate_vertices(data)
    _validate_edges(data)
    _validate_codon_set(data)
end


# validate codon_set contents
function _validate_codon_set(data::CodonGraphData)
    # do not allow empty codon sets
    length(data.codon_set) == 0 && throw(ArgumentError("\"codon_set\" is empty."))

    # do not allow duplicates
    length(data.codon_set) == length(Set(data.codon_set)) ||
        throw(ArgumentError("\"codon_set\" contains duplicate codons!"))

    # check if all codons in codon_set are valid
    for codon in data.codon_set
        codon in ALL_CODONS || throw(ArgumentError("codon \"$codon\" in \"codon_set\" is not a valid codon!"))
    end
end


# validate edge labels against edge count
function _validate_edges(data::CodonGraphData)
    # check if edge_labels is empty
    isempty(data.edge_labels) &&
        throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" is empty."))
    # check if graph contains edges
    ne(data.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no edges."))
    # check if edge_labels length matches number of edges in graph
    length(data.edge_labels) == ne(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"edge_labels\" length does not match number of edges in graph.",
        ),
    )

    # check if duplicates exist in edge_labels
    length(data.edge_labels) == length(Set(data.edge_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" contains duplicate labels."))

    # check if every label in edge_labels corresponds to an actual edge in the graph
    for (src_label, dst_label) in data.edge_labels
        # check if src_label exists in vert_idxs
        haskey(data.vert_idxs, src_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if dst_label exists in vert_idxs
        haskey(data.vert_idxs, dst_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if src_label exists in vert_labels
        src_label in data.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if dst_label exists in vert_labels
        dst_label in data.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if edge from src_label to dst_label exists in graph
        _has_edge_label(data, src_label, dst_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains edge from \"$src_label\" to \"$dst_label\" that does not exist in graph.",
            ),
        )
    end

    # check if every edge in graph has a corresponding label in edge_labels
    for edge in edges(data.graph)
        src_label = data.vert_labels[src(edge)]
        dst_label = data.vert_labels[dst(edge)]
        (src_label, dst_label) in data.edge_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: edge from \"$src_label\" to \"$dst_label\" does not have a corresponding label in \"edge_labels\".",
            ),
        )
    end
end


# validate vertices and index mappings
function _validate_vertices(data::CodonGraphData)
    # check if vert_labels is empty
    isempty(data.vert_labels) &&
        throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" is empty."))
    # check if vert_idxs is empty
    isempty(data.vert_idxs) && throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" is empty."))
    # check if graph contains vertices
    nv(data.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no vertices."))
    # check if no duplicates exist in vert_labels
    length(data.vert_labels) == length(Set(data.vert_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" contains duplicate labels."))

    # check if vert_labels only contains valid labels (length 1 or 2, only A,C,G,T or combinations thereof)
    for label in data.vert_labels
        (length(label) in (1, 2)) || throw(
            ArgumentError("inconsistent CodonGraphData: \"vert_labels\" contains invalid label \"$label\"."),
        )
        for char in label
            char in ALLOWED_BASES_STR || throw(
                ArgumentError(
                    "inconsistent CodonGraphData: \"vert_labels\" contains invalid label \"$label\".",
                ),
            )
        end
    end

    # check if vert_labels length matches graph size
    length(data.vert_labels) == nv(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_labels\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs length matches graph size
    length(data.vert_idxs) == nv(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_idxs\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs keys match vert_labels and indices are valid
    for (label, idx) in data.vert_idxs
        (1 <= idx <= nv(data.graph)) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"vert_idxs\" contains invalid vertex index \"$idx\".",
            ),
        )
        data.vert_labels[idx] == label ||
            throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" and vert_labels do not match."))
    end
end
