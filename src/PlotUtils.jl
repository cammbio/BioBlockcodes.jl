using CairoMakie
using GraphMakie
using NetworkLayout
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
    show_codon_graph(data::CodonGraphData; debug::Bool = false)

Create and display a plot representing `data`.

# Arguments

  - `data::CodonGraphData`: Graph data to visualize.

# Keyword Arguments

  - `debug::Bool`: Whether to emit debug logs.

# Returns

  - `Nothing`: Displays a figure.

# Throws

  - ArgumentError: If length of `all_vertex_labels` is not equal to number of vertices.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)
show_graph(data)
```
"""
function show_codon_graph(data::CodonGraphData; debug::Bool = false)
    # check if length of all_vertex_labels is equal to number of vertices
    length(data.all_vertex_labels) != nv(data.graph) && throw(
        ArgumentError(
            "Length of all_vertex_labels ($(length(data.all_vertex_labels))) is not equal to number of vertices ($(nv(data.graph))).",
        ),
    )

    debug && @debug "Showing graph..."
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
    ax.title = data.plot_title
    debug && @debug "all_vertex_labels in graph: $(data.all_vertex_labels)"
    graphplot!(
        ax,
        data.graph;
        layout = Spring(C = 50.0),
        nlabels = data.all_vertex_labels,
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
    display(fig)
end


"""
    show_multiple_codon_graphs(data_list::Vector{CodonGraphData}; debug::Bool = false)

Create and display a grid of plots for `data_list`.

# Arguments

  - `data_list::Vector{CodonGraphData}`: Graph data list to visualize.

# Keyword Arguments

  - `debug::Bool`: Whether to emit debug logs.

# Returns

  - `Nothing`: Displays a figure.

# Throws

  - ArgumentError: If `data_list` is empty.
  - ArgumentError: If length of `all_vertex_labels` is not equal to number of vertices for any data.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)
show_multiple_graphs([data])
```
"""
function show_multiple_codon_graphs(
    data_list::Vector{CodonGraphData};
    fig_title::Union{String, Nothing} = nothing,
    debug::Bool = false,
)
    # do not allow empty data_list
    isempty(data_list) && throw(ArgumentError("data_list cannot be empty."))
    # check if length of all_vertex_labels is equal to number of vertices for each data
    for data in data_list
        length(data.all_vertex_labels) != nv(data.graph) && throw(
            ArgumentError(
                "Length of all_vertex_labels ($(length(data.all_vertex_labels))) is not equal to number of vertices ($(nv(data.graph))).",
            ),
        )
    end

    debug && @debug "Showing multiple graphs..."
    # get grid size
    col_count = _get_col_count(length(data_list))
    # einheitliche Titelgröße anhand des längsten Titels schätzen
    max_title_len = maximum(length.(getfield.(data_list, :plot_title)))
    title_width_px = 1800 / col_count
    fontsize_est = title_width_px / max(max_title_len, 1) / 0.6
    uniform_titlesize = clamp(fontsize_est, 8, 26)

    # create plot figure with minimal padding
    fig = Figure(size = (1800, 900), figure_padding = 2)
    # set plot_title as label above all graphs if provided
    if !isempty(fig_title) && fig_title !== nothing
        Label(fig[0, 1:col_count], fig_title; fontsize = 24)
    end

    # add each graph to figure
    for (index, data) in enumerate(data_list)
        quotient, remainder = divrem(index - 1, col_count)
        row = quotient + 1
        column = remainder + 1
        debug && @debug "index: $index, row: $row, column: $column"
        ax = Axis(
            fig[row, column];
            xgridvisible = false,
            ygridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
        )
        ax.titlesize = uniform_titlesize
        ax.title = data.plot_title
        debug && @debug "all_vertex_labels in graph: $(data.all_vertex_labels)"
        graphplot!(
            ax,
            data.graph;
            layout = Spring(C = 50.0),
            nlabels = data.all_vertex_labels,
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
    # display(fig)

    return fig
end


# determine number of columns for grid layout based on amount of graphs
function _get_col_count(graph_count::Int)
    colum_count = max(1, ceil(Int, sqrt(graph_count)))
    return colum_count
end
