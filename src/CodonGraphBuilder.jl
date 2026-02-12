# DEPRECATED
using CairoMakie
using GraphMakie
using Graphs
function _add_vert_by_label!(data::CodonGraphData, label::String)
    # do not allow empty labels
    isempty(label) && throw(ArgumentError("label cannot be empty!"))
    # only allow labels of length 1 or 2
    !(length(label) in (1, 2)) &&
        throw(ArgumentError("label must be of length 1 or 2, got length $(length(label))!"))
    if label in data.vert_labels # vertex already exists
        return false
    else
        add_vertex!(data.graph)
        push!(data.vert_labels, label)
        data.vert_idxs[label] = nv(data.graph)
        return true
    end
end


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


# checks if a vertex with given label exists in the graph
function _has_vert_label(vert_idxs::Dict{String, Int}, label::String)
    return haskey(vert_idxs, label)
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