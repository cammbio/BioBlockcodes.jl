"""
    is_c3(cgd::CodonGraphData) -> Bool

Checks whether a codon graph satisfies the C3 property.

# Arguments

  - `cgd::CodonGraphData`: Codon graph data object that contains the graph.

# Returns

  - `Bool`: `true` if `cgd` is C3, otherwise `false`.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using GCATCodes

julia> codon_set = GCATCodes.LongDNA{4}.(["ATG", "TGA"]);

julia> cgd = CodonGraphData(codon_set);

julia> is_c3(cgd)
false
```
"""
function is_c3(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)

    # check if original graph is circular
    if is_circular(cgd)
        # check if shifted graphs are circular
        for shift in (1, 2)
            shifted_codon_set = left_shift_codon_set(cgd.codon_set, shift)
            shifted_cgd = CodonGraphData(shifted_codon_set)
            if !is_circular(shifted_cgd)
                return false
            end
        end
    else
        return false
    end
    return true
end


"""
    is_circular(cgd::CodonGraphData) -> Bool

Checks whether a codon graph is acyclic aka. contains no cycles.

# Arguments

  - `cgd::CodonGraphData`: Codon graph data object that contains the graph.

# Returns

  - `Bool`: `true` if no cycles are present, otherwise `false`.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using GCATCodes

julia> codon_set = GCATCodes.LongDNA{4}.(["AAC", "GTT"]);

julia> cgd = CodonGraphData(codon_set);

julia> is_circular(cgd)
true
```
"""
function is_circular(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)

    # check if graph has cycles longer than 0 aka. any cycles at all
    if _has_cycle_longer_than(cgd.graph, 0)
        return false
    end

    return true
end

function is_circular(codons::Vector{LongDNA{4}})
    cgd = CodonGraphData(codons)
    return is_circular(cgd)
end

"""
    is_comma_free(cgd::CodonGraphData) -> Bool

Checks whether a codon graph is comma-free.

# Arguments

  - `cgd::CodonGraphData`: Codon graph data object to check.

# Returns

  - `Bool`: `true` if `cgd` is comma-free, otherwise `false`.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using GCATCodes

julia> codon_set = GCATCodes.LongDNA{4}.(["CGA", "TAC"]);

julia> cgd = CodonGraphData(codon_set);

julia> is_comma_free(cgd)
true
```
"""
function is_comma_free(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)

    max_depth = 2
    # perform depth-limited DFS for each vertice
    for vertex in vertices(cgd.graph)
        if _dfs_depth(cgd.graph, vertex, 0, max_depth)
            return false
        end
    end

    return true
end


"""
    is_self_complementary(cgd::CodonGraphData) -> Bool

Checks whether a codon graph is invariant under complement and reversal.

# Arguments

  - `cgd::CodonGraphData`: Codon graph data object to check.

# Returns

  - `Bool`: `true` if the graph is self-complementary, otherwise `false`.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using GCATCodes

julia> codon_set = GCATCodes.LongDNA{4}.(["AGC", "CTG", "TGT", "ATC"]);

julia> cgd = CodonGraphData(codon_set);

julia> is_self_complementary(cgd)
false
```
"""
function is_self_complementary(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)

    # create complemented, reversed codon set
    comp_rev_codon_set = get_comp_rev_codon_set(cgd.codon_set)

    # create CodonGraphData from complemented, reversed codon set
    comp_rev_cgd = CodonGraphData(comp_rev_codon_set)
    _validate_cgd(comp_rev_cgd)

    # compare original graph with complemented, reversed graph
    if _is_codon_graphs_equal(cgd, comp_rev_cgd)
        return true
    else
        return false
    end
end


# check if a codon graph is strong C3 (i.e., C3 and expanded graph only has cycles of length 2)
"""
    is_strong_c3(cgd::CodonGraphData) -> Bool

Checks whether a codon graph is strong C3.

# Arguments

  - `cgd::CodonGraphData`: Codon graph data object to check.

# Returns

  - `Bool`: `true` if `cgd` is strong C3, otherwise `false`.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using GCATCodes

julia> codon_set = GCATCodes.LongDNA{4}.(["GGA", "TAA"]);

julia> cgd = CodonGraphData(codon_set);

julia> is_strong_c3(cgd)
true
```
"""
function is_strong_c3(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)

    # check if original graph is C3
    if is_c3(cgd)
        # create graph copy
        cgd_expanded = CodonGraphData(cgd.codon_set)
        _validate_cgd(cgd_expanded)
        # expand copy graph
        _expand_graph!(cgd_expanded)
        # validate expanded cgd object
        _validate_cgd(cgd_expanded)

        # check if expanded graph has cycles longer than 2
        max_length = 2
        if _has_cycle_longer_than(cgd_expanded.graph, max_length)
            return false
        end

        return true
    else
        return false
    end
end


# adds edge to graph by vertex labels
function _add_edge_by_label!(cgd::CodonGraphData, edge_label::Tuple{String, String})
    src_label, dst_label = edge_label
    # validate cgd
    _validate_cgd(cgd)
    # validate labels
    _validate_label(src_label)
    _validate_label(dst_label)
    # check if labels exist in vert_labels
    src_label in cgd.vert_labels ||
        throw(ArgumentError("vertex with label $src_label does not exist in vert_labels ."))
    dst_label in cgd.vert_labels ||
        throw(ArgumentError("vertex with label $dst_label does not exist in vert_labels."))
    # check if labels exist in vert_idxs
    haskey(cgd.vert_idxs, src_label) ||
        throw(ArgumentError("Vertex with label $src_label does not exist in vert_idxs."))
    haskey(cgd.vert_idxs, dst_label) ||
        throw(ArgumentError("Vertex with label $dst_label does not exist in vert_idxs."))

    # edge does not already exist
    if !_has_edge_label(cgd, edge_label)
        edge_idx = (cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])
        # add edge to graph and update affected cgd fields
        if add_edge!(cgd.graph, edge_idx)
            push!(cgd.edge_labels, edge_label)
            return true
        else
            return false
        end
    else
        return false
    end
end


# adds vertex to graph by label
function _add_vert_by_label!(cgd::CodonGraphData, vert_label::String)
    # validate cgd
    _validate_cgd(cgd)
    # validate label
    _validate_label(vert_label)

    # vertex does not already exist
    if !_has_vert_label(cgd, vert_label)
        if add_vertex!(cgd.graph)
            push!(cgd.vert_labels, vert_label)
            cgd.vert_idxs[vert_label] = nv(cgd.graph)
            return true
        else
            return false
        end
    else
        return false
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


# expand graph by adding N2 and N3N1 for each codon N1N2N3 as vertices and connect them accordingly (if possible)
function _expand_graph!(cgd::CodonGraphData)
    # validate cgd
    _validate_cgd(cgd)
    for codon in cgd.codon_set
        n2 = string(codon[2])
        n3n1 = string(codon[3], codon[1])
        _add_vert_by_label!(cgd, n2)
        _add_vert_by_label!(cgd, n3n1)
        _add_edge_by_label!(cgd, (n2, n3n1))
        _add_edge_by_label!(cgd, (n3n1, n2))
    end
end


# check if graph has cycle longer than max_length
function _has_cycle_longer_than(graph::SimpleDiGraph, max_length::Int)
    # validate graph
    _validate_graph(graph)
    # do not allow negative max_length
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
function _has_edge_label(cgd::CodonGraphData, edge_label::Tuple{String, String})
    src_label, dst_label = edge_label
    # validate labels
    _validate_label(src_label)
    _validate_label(dst_label)
    # check if labels exist in vert_labels
    src_label in cgd.vert_labels ||
        throw(ArgumentError("vertex with label $src_label does not exist in vert_labels."))
    dst_label in cgd.vert_labels ||
        throw(ArgumentError("vertex with label $dst_label does not exist in vert_labels."))
    # check if labels exist in vert_idxs
    haskey(cgd.vert_idxs, src_label) ||
        throw(ArgumentError("vertex with label $src_label does not exist in vert_idxs."))
    haskey(cgd.vert_idxs, dst_label) ||
        throw(ArgumentError("vertex with label $dst_label does not exist in vert_idxs."))

    src_index = cgd.vert_idxs[src_label]
    dst_index = cgd.vert_idxs[dst_label]
    return (edge_label in cgd.edge_labels) && has_edge(cgd.graph, src_index, dst_index)
end


# checks if a vertex with given label exists in the graph
function _has_vert_label(cgd::CodonGraphData, label::String)
    # validate label
    _validate_label(label)

    return label in cgd.vert_labels &&
           haskey(cgd.vert_idxs, label) &&
           has_vertex(cgd.graph, cgd.vert_idxs[label])
end


# check if two codon graphs are equal by comparing vertices and edges
function _is_codon_graphs_equal(cgd_1::CodonGraphData, cgd_2::CodonGraphData)
    # check if same amount of vertices
    if nv(cgd_1.graph) != nv(cgd_2.graph)
        return false
    end
    # check if same amount of edges
    if ne(cgd_1.graph) != ne(cgd_2.graph)
        return false
    end

    # check if same vertice labels
    for index in 1:nv(cgd_1.graph)
        if !_has_vert_label(cgd_2, cgd_1.vert_labels[index])
            return false
        end
    end

    # check if same edges
    for edge in edges(cgd_1.graph)
        edge_label = (cgd_1.vert_labels[src(edge)], cgd_1.vert_labels[dst(edge)])
        if _has_edge_label(cgd_2, edge_label)
            continue
        else
            return false
        end
    end

    return true
end



