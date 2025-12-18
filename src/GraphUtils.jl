using CairoMakie
using GraphMakie
using Graphs
# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# construct graph from data
function construct_graph!(data::CodonGraphData; show_debug::Bool = false, show_plot::Bool = false)
    # check if any duplicates in codons
    if length(data.codon_set) == length(unique(data.codon_set)) # no duplicates if true
        show_debug && @debug """
        Before adding vertices and edges:
        graph: $(data.graph)
        codon_set: $(data.codon_set)
        all_vertex_labels: $(data.all_vertex_labels)
        base_vertex_labels: $(data.base_vertex_labels)
        added_vertex_labels: $(data.added_vertex_labels)
        all_edge_labels: $(data.all_edge_labels)
        base_edge_labels: $(data.base_edge_labels)
        added_edge_labels: $(data.added_edge_labels)
        vertex_index: $(data.vertex_index)
        """

        # extract vertice labels from codon set and add vertices to graph
        create_vertices!(data)
        # create mapping from vertice label to vertice index in graph
        data.vertex_index = Dict(label => index for (index, label) in enumerate(data.base_vertex_labels))
        # connect edges
        connect_edges!(data)

        show_debug && @debug """
        After adding vertices and edges:
        graph: $(data.graph)
        codon_set: $(data.codon_set)
        all_vertex_labels: $(data.all_vertex_labels)
        base_vertex_labels: $(data.base_vertex_labels)
        added_vertex_labels: $(data.added_vertex_labels)
        all_edge_labels: $(data.all_edge_labels)
        base_edge_labels: $(data.base_edge_labels)
        added_edge_labels: $(data.added_edge_labels)
        vertex_index: $(data.vertex_index)
        """
    end
    show_debug && @debug "Graph construction from codon set finished: $(data.codon_set)"
    if show_plot
        show_graph(data; show_debug = show_debug)
    end
end


# create needed vertices for the graph by iterating through codons and collecting all singular and tuple bases
function create_vertices!(data::CodonGraphData; show_debug::Bool = false)
    # use a temporary set to avoid duplicates and increase lookup speed
    temp_labels = Set{String}()

    # iterate through codon set and extract needed vertice labels
    for codon in data.codon_set
        # get first and last character of codon
        push!(temp_labels, string(codon[1])) # first base
        push!(temp_labels, string(codon[3])) # third base
        # get first two and last two characters of codon
        push!(temp_labels, string(codon[1:2])) # first tuple
        push!(temp_labels, string(codon[2:3])) # second tuple
    end

    # sort and copy temp_labels to all_vertex_labels and base_vertex_labels fields
    labels = collect(temp_labels)
    sort!(labels, by = x -> (length(x), x))
    data.all_vertex_labels = labels
    data.base_vertex_labels = copy(labels)

    # add vertices to graph
    for _ in 1:length(data.all_vertex_labels)
        add_vertex!(data.graph)
    end
end


# connect the vertices in the graph based on the codons by connecting the first base to the second tuple and
# the first tuple to the third base
function connect_edges!(data::CodonGraphData; show_debug::Bool = false)
    graph = data.graph
    vertex_index = data.vertex_index
    base_vertex_labels = data.base_vertex_labels
    all_edge_labels = data.all_edge_labels
    base_edge_labels = data.base_edge_labels

    # iterate through codon set and add edges to graph
    for codon in data.codon_set
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