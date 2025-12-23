using CairoMakie
using GraphMakie
using NetworkLayout
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
    show_graph(data::CodonGraphData; show_debug::Bool = false)

Create and display a plot representing `data`.

# Arguments
- `data::CodonGraphData`: Graph data to visualize.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Nothing`: Displays a figure.

# Throws
- None.

# Example
```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph!(data)
show_graph(data)
```
"""
function show_graph(data::CodonGraphData; show_debug::Bool = false)
    show_debug && @debug "Showing graph..."
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
    # combine vertice labels with manually added vertice labels
    show_debug && @debug "all_vertex_labels in graph: $(data.all_vertex_labels)"
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
        node_size = 50,
        arrow_shift = :end,
        arrow_size = 12,
        edge_width = 2,
        edge_curvature = 0.9,
    )
    display(fig)
end


"""
    show_multiple_graphs(data_list::Vector{CodonGraphData}; show_debug::Bool = false)

Create and display a grid of plots for `data_list`.

# Arguments
- `data_list::Vector{CodonGraphData}`: Graph data list to visualize.

# Keyword Arguments
- `show_debug::Bool`: Whether to emit debug logs.

# Returns
- `Nothing`: Displays a figure.

# Throws
- None.

# Example
```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph!(data)
show_multiple_graphs([data])
```
"""
function show_multiple_graphs(data_list::Vector{CodonGraphData}; show_debug::Bool = false)
    show_debug && @debug "Showing multiple graphs..."
    # get grid size
    amount_graphs = length(data_list)
    number_rows, number_columns = _grid_size(amount_graphs)
    # create plot figure
    fig = Figure(size = (1800, 900))
    # add each graph to figure
    for (index, data) in enumerate(data_list)
        row = fld(index - 1, number_columns) + 1
        column = mod(index - 1, number_columns) + 1
        show_debug && @debug "index: $index, row: $row, column: $column"
        ax = Axis(
            fig[row, column];
            xgridvisible = false,
            ygridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
        )
        ax.title = data.plot_title
        # combine vertice labels with manually added vertice labels
        combine_vertice_labels = vcat(data.vertice_labels, data.added_vertice_labels)
        show_debug && @debug "all_vertex_labels in graph: $(data.all_vertex_labels)"
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
            node_size = 50,
            arrow_shift = :end,
            arrow_size = 12,
            edge_width = 2,
            edge_curvature = 0.9,
        )
    end
    display(fig)
end


"""
    _grid_size(amount_graphs::Int) -> Tuple{Int, Int}

Return the grid size (number_rows, number_columns) for `amount_graphs`.

# Arguments
- `amount_graphs::Int`: Number of graphs to place.

# Keyword Arguments
- None.

# Returns
- `Tuple{Int, Int}`: `(number_rows, number_columns)` grid size.

# Throws
- None.

# Example
```julia
_grid_size(5)
```
"""
function _grid_size(amount_graphs::Int)
    number_columns = max(1, ceil(Int, sqrt(amount_graphs)))
    number_rows = ceil(Int, amount_graphs / number_columns)
    return number_rows, number_columns
end
