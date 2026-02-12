using BioSequences
using Graphs


# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    graph::Graphs.SimpleDiGraph
    codon_set::Vector{LongDNA{4}}
    vert_labels::Vector{String}
    edge_labels::Vector{Tuple{String, String}}
    vert_idxs::Dict{String, Int}
    graph_title::String
end


# constructor for CodonGraphData that takes a codon set and builds the other fields
function CodonGraphData(codon_set::Vector{LongDNA{4}}; graph_title::String = "")
    isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty!"))
    any(length(codon) != 3 for codon in codon_set) &&
        throw(ArgumentError("All codons in codon set must have length 3!"))
    length(codon_set) != length(Set(codon_set)) &&
        throw(ArgumentError("Codon set cannot contain duplicate codons!"))
    valid_bases = Set((DNA_A, DNA_C, DNA_G, DNA_T))
    any(any(!(base in valid_bases) for base in codon) for codon in codon_set) &&
        throw(ArgumentError("Codon set contains invalid DNA bases! Only A, C, G, T are allowed."))

    obj = CodonGraphData(
        Graphs.SimpleDiGraph(0),
        copy(codon_set),
        String[],
        Tuple{String, String}[],
        Dict{String, Int}(),
        graph_title,
    )

    # build graph by adding vertices and edges based on codon set
    _add_vertices!(obj)
    obj.vert_idxs = Dict(label => index for (index, label) in enumerate(obj.vert_labels))
    _add_edges!(obj)

    return obj
end


# adds edges to graph after extracting needed labels from codon set
function _add_edges!(obj::CodonGraphData)
    # iterate through codon set and add edges to graph
    for codon in obj.codon_set
        # get needed vertex IDs
        first_base_idx = obj.vert_idxs[string(codon[1])]
        third_base_idx = obj.vert_idxs[string(codon[3])]
        first_tuple_idx = obj.vert_idxs[string(codon[1:2])]
        second_tuple_idx = obj.vert_idxs[string(codon[2:3])]

        # add edge_labels to edge_labels fields
        push!(obj.edge_labels, (obj.vert_labels[first_base_idx], obj.vert_labels[second_tuple_idx]))
        push!(obj.edge_labels, (obj.vert_labels[first_tuple_idx], obj.vert_labels[third_base_idx]))

        # add edges to graph
        add_edge!(obj.graph, first_base_idx, second_tuple_idx)
        add_edge!(obj.graph, first_tuple_idx, third_base_idx)
    end
end


# adds vertices to graph after extracting needed labels from codon set
function _add_vertices!(obj::CodonGraphData)
    # use a temporary set to avoid duplicates and increase lookup speed
    temp_labels = Set{String}()

    # iterate through codon set and extract needed vertice labels
    for codon in obj.codon_set
        # get first and last character of codon
        push!(temp_labels, string(codon[1]))
        push!(temp_labels, string(codon[3]))
        # get first and second tuple of codon
        push!(temp_labels, string(codon[1:2]))
        push!(temp_labels, string(codon[2:3]))
    end

    # sort and copy temp_labels to vert_labels
    labels::Vector{String} = collect(temp_labels)
    sort!(labels, by = x -> (length(x), x))
    append!(obj.vert_labels, labels)

    # add a vertex for each label to graph
    for _ in 1:length(obj.vert_labels)
        add_vertex!(obj.graph)
    end
end
