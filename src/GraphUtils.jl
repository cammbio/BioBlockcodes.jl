using CairoMakie
using GraphMakie
using Graphs
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
    construct_graph_data!(data::CodonGraphData; show_debug::Bool = false) -> Bool

Construct the graph for `data` from its codon set.

# Arguments

  - `data::CodonGraphData`: Graph data to populate.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `True`: If graph construction was successful.

# Throws

  - `ArgumentError`: If a `data` field is not empty before construction (except codon_set) or if the graph is not empty.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data; show_debug = true)
```
"""
function construct_graph_data!(data::CodonGraphData; show_debug::Bool = false)
    show_debug &&
        @debug "Debug logs for construct_graph_data!-------------------------------------------------------------------------------------------------------------------------------------------"
    # do not allow populated data fields (only codon_set is allowed to be populated)
    _check_data_fields_empty(data, show_debug = show_debug)

    show_debug && @debug """Before adding vertices and edges:
    graph: $(data.graph)
    codon_set: $(data.codon_set)
    all_vertex_labels: $(data.all_vertex_labels)
    base_vertex_labels: $(data.base_vertex_labels)
    added_vertex_labels: $(data.added_vertex_labels)
    all_edge_labels: $(data.all_edge_labels)
    base_edge_labels: $(data.base_edge_labels)
    added_edge_labels: $(data.added_edge_labels)
    vertex_index: $(data.vertex_index)"""

    # add vertices to graph based on codon set
    _add_vertices_by_codon_set!(
        data.graph,
        data.codon_set,
        data.all_vertex_labels,
        data.base_vertex_labels;
        show_debug = show_debug,
    )

    # create mapping from vertice label to vertice index in graph
    data.vertex_index = Dict(label => index for (index, label) in enumerate(data.base_vertex_labels))

    # add edges to graph based on codon set
    _add_edges_by_codon_set!(
        data.graph,
        data.codon_set,
        data.base_vertex_labels,
        data.all_edge_labels,
        data.base_edge_labels,
        data.vertex_index;
        show_debug = show_debug,
    )

    show_debug && @debug """After adding vertices and edges:
    graph: $(data.graph)
    codon_set: $(data.codon_set)
    all_vertex_labels: $(data.all_vertex_labels)
    base_vertex_labels: $(data.base_vertex_labels)
    added_vertex_labels: $(data.added_vertex_labels)
    all_edge_labels: $(data.all_edge_labels)
    base_edge_labels: $(data.base_edge_labels)
    added_edge_labels: $(data.added_edge_labels)
    vertex_index: $(data.vertex_index)
    Graph construction from codon set successfully finished: $(data.codon_set)"""
    return true
end


"""
    add_vertex_by_label!(data::CodonGraphData, label::String; show_debug::Bool = false) -> Bool

Add a vertex with `label` to the graph if it does not already exist.

# Arguments

  - `data::CodonGraphData`: Graph data to modify.
  - `label::String`: Vertex label to add.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if the vertex was added, otherwise `false`.

# Throws

  - `ArgumentError`: If `label` is empty.
  - `ArgumentError`: If `label` is not of length 1 or 2.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)
add_vertex_by_label!(data, "AG")
```
"""
function add_vertex_by_label!(data::CodonGraphData, label::String; show_debug::Bool = false)
    show_debug &&
        @debug "Debug logs for add_vertex_by_label!-------------------------------------------------------------------------------------------------------------------------------------------"
    # do not allow empty labels
    isempty(label) && throw(ArgumentError("label cannot be empty!"))
    # only allow labels of length 1 or 2
    !(length(label) in (1, 2)) &&
        throw(ArgumentError("label must be of length 1 or 2, got length $(length(label))!"))
    if label in data.all_vertex_labels # vertex already exists
        show_debug && @debug "Vertice $label already exists in graph -> not added."
        return false
    else # vertex does not already exist
        # update affected data fields
        add_vertex!(data.graph) # add vertex to graph
        push!(data.added_vertex_labels, label) # add to manually added vertex labels
        push!(data.all_vertex_labels, label) # add to all vertex labels
        data.vertex_index[label] = nv(data.graph) # map label to vertex index
        show_debug && @debug "Added vertex: $label"
        return true
    end
end


"""
    add_edge_by_label!(
        data::CodonGraphData,
        src_label::String,
        dst_label::String;
        show_debug::Bool = false,
    ) -> Bool

Add a labeled edge if it does not already exist.

# Arguments

  - `data::CodonGraphData`: Graph data to modify.
  - `from_label::String`: Source label.
  - `to_label::String`: Destination label.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Bool`: `true` if the edge was added, otherwise `false`.

# Throws

  - `ArgumentError`: If `src_label` or `dst_label` is empty.
  - `ArgumentError`: If `src_label` or `dst_label` does not exist in the graph.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)
add_edge_by_label!(data, "C", "GT")
```
"""
function add_edge_by_label!(
    data::CodonGraphData,
    src_label::String,
    dst_label::String;
    show_debug::Bool = false,
)
    # do not allow empty labels or not existing labels in vertex_index
    isempty(src_label) && throw(ArgumentError("from_label cannot be empty!"))
    isempty(dst_label) && throw(ArgumentError("to_label cannot be empty!"))
    !haskey(data.vertex_index, src_label) &&
        throw(ArgumentError("Vertex with label $src_label does not exist in graph!"))
    !haskey(data.vertex_index, dst_label) &&
        throw(ArgumentError("Vertex with label $dst_label does not exist in graph!"))

    if _has_edge_label(data, src_label, dst_label; show_debug = show_debug) # edge already exists
        show_debug && @debug "Edge $src_label -> $dst_label already exists in graph -> not added."
        return false
    else # edge does not already exist
        # get vertex indices from labels
        from_index = data.vertex_index[src_label]
        to_index = data.vertex_index[dst_label]

        # add edge to graph and update affected data fields
        push!(data.added_edge_labels, (src_label, dst_label))
        push!(data.all_edge_labels, (src_label, dst_label))
        add_edge!(data.graph, from_index, to_index)
        show_debug && @debug "Added edge: $src_label -> $dst_label"
        return true
    end
end


# ---------------------------------------------- HELPERS ----------------------------------------------
# adds vertices to graph after extracting needed labels from codon set
function _add_vertices_by_codon_set!(
    graph::Graphs.SimpleDiGraph,
    codon_set::Vector{LongDNA{4}},
    all_vertex_labels::Vector{String},
    base_vertex_labels::Vector{String};
    show_debug::Bool = false,
)
    # use a temporary set to avoid duplicates and increase lookup speed
    temp_labels = Set{String}()

    # iterate through codon set and extract needed vertice labels
    for codon in codon_set
        # get first and last character of codon
        push!(temp_labels, string(codon[1])) # first base
        push!(temp_labels, string(codon[3])) # third base
        # get first and second tuple of codon
        push!(temp_labels, string(codon[1:2])) # first tuple
        push!(temp_labels, string(codon[2:3])) # second tuple
    end

    # sort and copy temp_labels to all_vertex_labels and base_vertex_labels
    labels::Vector{String} = collect(temp_labels) # turn set into vector
    sort!(labels, by = x -> (length(x), x)) # sort by length then lexicographically
    append!(all_vertex_labels, labels)
    append!(base_vertex_labels, labels)

    # add a vertex for each label to graph
    for _ in 1:length(base_vertex_labels)
        add_vertex!(graph)
    end
end


# adds edges to graph after extracting needed labels from codon set
function _add_edges_by_codon_set!(
    graph::Graphs.SimpleDiGraph,
    codon_set::Vector{LongDNA{4}},
    base_vertex_labels::Vector{String},
    all_edge_labels::Vector{Tuple{String, String}},
    base_edge_labels::Vector{Tuple{String, String}},
    vertex_index::Dict{String, Int};
    show_debug::Bool = false,
)
    # iterate through codon set and add edges to graph
    for codon in codon_set
        # add_edge_by_codon(graph, codon; show_debug = show_debug)
        # get needed vertex IDs
        first_base_id = vertex_index[string(codon[1])]
        third_base_id = vertex_index[string(codon[3])]
        first_tuple_id = vertex_index[string(codon[1:2])]
        second_tuple_id = vertex_index[string(codon[2:3])]

        # add edge_labels to all_edge_labels and base_edge_labels fields
        push!(all_edge_labels, (base_vertex_labels[first_base_id], base_vertex_labels[second_tuple_id]))
        push!(all_edge_labels, (base_vertex_labels[first_tuple_id], base_vertex_labels[third_base_id]))
        push!(base_edge_labels, (base_vertex_labels[first_base_id], base_vertex_labels[second_tuple_id]))
        push!(base_edge_labels, (base_vertex_labels[first_tuple_id], base_vertex_labels[third_base_id]))

        # add edges to graph
        add_edge!(graph, first_base_id, second_tuple_id)
        add_edge!(graph, first_tuple_id, third_base_id)
    end
end


# checks if a vertex with given label exists in the graph
function _has_vertex_label(vertex_index::Dict{String, Int}, label::String; show_debug::Bool = false)
    show_debug && @debug "Checking if vertice label $label exists in graph..."
    return haskey(vertex_index, label)
end


# checks if an edge with given labels exists in the graph
function _has_edge_label(data::CodonGraphData, src_label::String, dst_label::String; show_debug::Bool = false)
    # check if labels exist
    haskey(data.vertex_index, src_label) || return false
    haskey(data.vertex_index, dst_label) || return false

    src_index = data.vertex_index[src_label]
    dst_index = data.vertex_index[dst_label]

    return has_edge(data.graph, src_index, dst_index)
end


# checks if all data fields (except codon_set) are empty
function _check_data_fields_empty(data::CodonGraphData; show_debug::Bool = false)
    nv(data.graph) != 0 && throw(ArgumentError("Graph contains vertices."))
    ne(data.graph) != 0 && throw(ArgumentError("Graph contains edges."))
    !isempty(data.all_vertex_labels) && throw(ArgumentError("all_vertex_labels is not empty."))
    !isempty(data.base_vertex_labels) && throw(ArgumentError("base_vertex_labels is not empty."))
    !isempty(data.added_vertex_labels) && throw(ArgumentError("added_vertex_labels is not empty."))
    !isempty(data.all_edge_labels) && throw(ArgumentError("all_edge_labels is not empty."))
    !isempty(data.base_edge_labels) && throw(ArgumentError("base_edge_labels is not empty."))
    !isempty(data.added_edge_labels) && throw(ArgumentError("added_edge_labels is not empty."))
    !isempty(data.vertex_index) && throw(ArgumentError("vertex_index is not empty."))
    return true
end