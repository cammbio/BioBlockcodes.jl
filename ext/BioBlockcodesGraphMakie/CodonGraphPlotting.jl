# determine number of columns for grid layout based on amount of graphs
function _get_col_count(graph_count::Int)
    # do not allow 0 or negative graph count
    if graph_count <= 0
        throw(ArgumentError("graph_count must be a positive integer."))
    end

    col_count = max(1, ceil(Int, sqrt(graph_count)))
    return col_count
end


"""
    plot_codon_graph(cgd::CodonGraphData; fig_size::Tuple{Int, Int} = (1800, 900)) -> Figure

Creates a single visualization of a codon graph in a Figure.

# Arguments

  - `cgd::CodonGraphData`: Codon graph to plot.

# Keyword Arguments

  - `fig_size::Tuple{Int, Int}=(1800, 900)`: Size of the generated figure in pixels.

# Returns

  - `Figure`: Makie figure containing the plotted graph.

# Throws

  - `ArgumentError`: If `cgd` is invalid.

# Examples

```jldoctest
julia> using BioBlockcodes, BioSequences

julia> codon_set = [dna"CAA", dna"GTC"];

julia> cgd = CodonGraphData(codon_set);

julia> fig = plot_codon_graph(cgd);

julia> typeof(fig)
Makie.Figure
```
"""
function plot_codon_graph(cgd::CodonGraphData; fig_size::Tuple{Int, Int} = (600, 400), spring = 100)
    # create plot figure
    fig = Figure(size = fig_size)
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
    ax.title = cgd.graph_title
    graphplot!(
        ax,
        cgd.graph;
        layout = NetworkLayout.Spring(C = spring),
        nlabels = cgd.vert_labels,
        nlabels_color = :black,
        nlabels_size = 28,
        nlabels_offset = Point2f(0, 0),
        nlabels_align = (:center, :center),
        node_color = :lightgrey,
        node_size = 50,
        arrow_shift = :end,
        arrow_size = 10,
        edge_width = 2,
        edge_curvature = 0.9,
    )

    return fig
end


"""
    plot_multiple_codon_graphs(cgd_list::Vector{CodonGraphData}; fig_title::Union{String, Nothing}=nothing, fig_size::Tuple{Int, Int}=(1800, 900)) -> Figure

Creates a grid visualization of multiple codon graphs in a Figure.

# Arguments

  - `cgd_list::Vector{CodonGraphData}`: List of codon graphs to plot.

# Keyword Arguments

  - `fig_title::Union{String, Nothing}=nothing`: Optional global title above the grid.
  - `fig_size::Tuple{Int, Int}=(1800, 900)`: Size of the generated figure in pixels.

# Returns

  - `Figure`: Makie figure containing all plotted graphs.

# Throws

  - `ArgumentError`: If `cgd_list` is empty or contains invalid entries.

# Examples

```jldoctest
julia> using BioBlockcodes, BioSequences

julia> codon_set_1 = [dna"ATA", dna"GGA"];

julia> codon_set_2 = [dna"TTC", dna"GCA"];

julia> cgd_1 = CodonGraphData(codon_set_1);

julia> cgd_2 = CodonGraphData(codon_set_2);

julia> cgd_list = [cgd_1, cgd_2];

julia> fig = plot_multiple_codon_graphs(cgd_list);

julia> typeof(fig)
Makie.Figure
```
"""
function plot_multiple_codon_graphs(
    cgd_list::Vector{CodonGraphData};
    fig_title::Union{String, Nothing} = nothing,
    fig_size::Tuple{Int, Int} = (1800, 900),
)
    # do not allow empty cgd_list
    isempty(cgd_list) && throw(ArgumentError("cgd_list is empty."))
    # validate cgd objects
    BioBlockcodes._validate_cgd.(cgd_list)

    # get number of columns for grid layout based on amount of graphs
    col_count = _get_col_count(length(cgd_list))
    # einheitliche Titelgröße anhand des längsten Titels schätzen
    max_title_len = maximum(length.(getfield.(cgd_list, :graph_title)))
    title_width = 1800 / col_count
    fontsize_est = title_width / max(max_title_len, 1) / 0.6
    title_size = clamp(fontsize_est, 8, 26)

    # create plot figure with minimal padding
    fig = Figure(size = fig_size, figure_padding = 2)
    # set fig_title as label above all graphs if provided
    if fig_title !== nothing && !isempty(fig_title)
        Label(fig[0, 1:col_count], fig_title; fontsize = 24)
    end

    # add each graph to figure
    for (idx, cgd) in enumerate(cgd_list)
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
        ax.title = cgd.graph_title
        graphplot!(
            ax,
            cgd.graph;
            layout = NetworkLayout.Spring(C = 50.0),
            nlabels = cgd.vert_labels,
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
