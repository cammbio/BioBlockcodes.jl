# using BioSequences
# using Graphs

# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    codon_set::Vector{LongDNA{4}}
    edge_labels::Vector{Tuple{String, String}}
    graph::Graphs.SimpleDiGraph
    graph_title::String
    vert_idxs::Dict{String, Int}
    vert_labels::Vector{String}
end


# constructor for CodonGraphData that takes a codon set and builds the other fields
"""
    CodonGraphData(codon_set::Vector{LongDNA{4}}; graph_title::String = "") -> CodonGraphData

Constructs a `CodonGraphData` instance from a codon set and builds the corresponding graph.

# Arguments
- `codon_set::Vector{LongDNA{4}}`: Input codons used to build the graph.

# Keyword Arguments
- `graph_title::String=""`: Optional title for the graph.

# Returns
- `CodonGraphData`: Fully initialized structure with vertices, edges, and labels.

# Throws
- `ArgumentError`: If `codon_set` is invalid or graph construction fails.

# Examples
```jldoctest
julia> using GCATCodes

julia> using BioSequences: LongDNA

julia> codon_set = LongDNA{4}.(["CCA", "GAT"]);

julia> cgd = CodonGraphData(codon_set);

julia> cgd isa CodonGraphData
true
```
"""
function CodonGraphData(codon_set::Vector{LongDNA{4}}; graph_title::String = "")
    # validate codon_set
    _validate_codon_set(codon_set)

    cgd = CodonGraphData(
        copy(codon_set),
        Tuple{String, String}[],
        Graphs.SimpleDiGraph(0),
        graph_title,
        Dict{String, Int}(),
        String[],
    )

    # build graph by adding vertices and edges based on codon set
    _add_vertices!(cgd)
    cgd.vert_idxs = Dict(label => index for (index, label) in enumerate(cgd.vert_labels))
    _add_edges!(cgd)

    return cgd
end


# adds edges to graph after extracting needed labels from codon set
function _add_edges!(cgd::CodonGraphData)
    # iterate through codon set and add edges to graph
    for codon in cgd.codon_set
        # get needed vertex IDs
        first_base_idx = cgd.vert_idxs[string(codon[1])]
        third_base_idx = cgd.vert_idxs[string(codon[3])]
        first_tuple_idx = cgd.vert_idxs[string(codon[1:2])]
        second_tuple_idx = cgd.vert_idxs[string(codon[2:3])]

        first_edge = (first_base_idx, second_tuple_idx)
        second_edge = (first_tuple_idx, third_base_idx)
        first_edge_label = (cgd.vert_labels[first_base_idx], cgd.vert_labels[second_tuple_idx])
        second_edge_label = (cgd.vert_labels[first_tuple_idx], cgd.vert_labels[third_base_idx])

        # add edges to graph and update fields
        if add_edge!(cgd.graph, first_edge[1], first_edge[2])
            push!(cgd.edge_labels, first_edge_label)
        else
            throw(
                ArgumentError(
                    "Failed to add edge from $(first_edge_label[1]) to $(first_edge_label[2]) for codon $codon. This should not happen if the graph is built correctly, so there may be an issue with the codon set or graph construction.",
                ),
            )
        end
        if add_edge!(cgd.graph, second_edge[1], second_edge[2])
            push!(cgd.edge_labels, second_edge_label)
        else
            throw(
                ArgumentError(
                    "Failed to add edge from $(second_edge_label[1]) to $(second_edge_label[2]) for codon $codon. This should not happen if the graph is built correctly, so there may be an issue with the codon set or graph construction.",
                ),
            )
        end
    end
end


# adds vertices to graph after extracting needed labels from codon set
function _add_vertices!(cgd::CodonGraphData)
    # use a temporary set to avoid duplicates
    temp_labels = Set{String}()

    # iterate through codon set and extract needed vertice labels
    for codon in cgd.codon_set
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
    append!(cgd.vert_labels, labels)

    # add a vertex for each label to graph
    for _ in 1:length(cgd.vert_labels)
        add_vertex!(cgd.graph)
    end
end



