using CairoMakie
using GraphMakie
using NetworkLayout


function plot_codon_graph(data::CodonGraphData)
    # check if length of vert_labels is equal to number of vertices
    length(data.vert_labels) != nv(data.graph) && throw(
        ArgumentError(
            "Length of vert_labels ($(length(data.vert_labels))) is not equal to number of vertices ($(nv(data.graph))).",
        ),
    )

    # create plot figure
    fig = Figure(size = (1800, 900))
    ax = Axis(
        fig[1, 1];
        xgridvisible = false,
        ygridvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
    )
    hidespines!(ax) # remove axis spines
    ax.title = data.graph_title
    graphplot!(
        ax,
        data.graph;
        layout = Spring(C = 50.0),
        nlabels = data.vert_labels,
        nlabels_color = :white,
        nlabels_size = 18,
        nlabels_offset = Point2f(0, 0),
        nlabels_align = (:center, :center),
        node_color = :black,
        node_size = 30,
        arrow_shift = :end,
        arrow_size = 12,
        edge_width = 2,
        edge_curvature = 0.9,
    )

    return fig
end


function plot_multiple_codon_graphs(
    data_list::Vector{CodonGraphData};
    fig_title::Union{String, Nothing} = nothing,
)
    # do not allow empty data_list
    isempty(data_list) && throw(ArgumentError("data_list cannot be empty."))
    # check if length of vert_labels is equal to number of vertices for each data
    for data in data_list
        length(data.vert_labels) != nv(data.graph) && throw(
            ArgumentError(
                "Length of vert_labels ($(length(data.vert_labels))) is not equal to number of vertices ($(nv(data.graph))).",
            ),
        )
    end

    # get grid size
    col_count = _get_col_count(length(data_list))
    # einheitliche Titelgröße anhand des längsten Titels schätzen
    max_title_len = maximum(length.(getfield.(data_list, :graph_title)))
    title_width = 1800 / col_count
    fontsize_est = title_width / max(max_title_len, 1) / 0.6
    title_size = clamp(fontsize_est, 8, 26)

    # create plot figure with minimal padding
    fig = Figure(size = (1800, 900), figure_padding = 2)
    # set fig_title as label above all graphs if provided
    if !isempty(fig_title) && fig_title !== nothing
        Label(fig[0, 1:col_count], fig_title; fontsize = 24)
    end

    # add each graph to figure
    for (idx, data) in enumerate(data_list)
        quot, rem = divrem(idx - 1, col_count)
        row = quot + 1
        column = rem + 1
        ax = Axis(
            fig[row, column];
            xgridvisible = false,
            ygridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
        )
        ax.titlesize = title_size
        ax.title = data.graph_title
        graphplot!(
            ax,
            data.graph;
            layout = Spring(C = 50.0),
            nlabels = data.vert_labels,
            nlabels_color = :white,
            nlabels_size = 18,
            nlabels_offset = Point2f(0, 0),
            nlabels_align = (:center, :center),
            node_color = :black,
            node_size = 30,
            arrow_shift = :end,
            arrow_size = 12,
            edge_width = 2,
            edge_curvature = 0.9,
        )
    end

    return fig
end


# determine number of columns for grid layout based on amount of graphs
function _get_col_count(graph_count::Int)
    col_count = max(1, ceil(Int, sqrt(graph_count)))
    return col_count
end
