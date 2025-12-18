module GCATCodes

# load other source files
include("Types.jl")
include("CyclesAnalysis.jl") # needs Types.jl
include("CodonUtils.jl") # needs Types.jl
include("PlotUtils.jl") # needs Types.jl
include("CodonGraphAnalysis.jl") # needs Types.jl and CodonUtils.jl
include("GraphUtils.jl") # needs Types.jl and PlotUtils.jl

# export types
export
# Types.jl
    CodonGraphData
#

# export functions
export
    # CodonGraphAnalysis.jl
    is_circular,
    is_comma_free,
    is_self_complementary,
    is_c3,
    is_graphs_identical,
    has_vertex_label,
    has_edge_label,
    # CodonUtils.jl
    get_complemented_reversed_codon_set,
    get_complemented_codon_set,
    get_complemented_codon,
    get_reversed_codon_set,
    get_complemented_codon,
    get_reversed_codon,
    get_complemented_base,
    left_shift_codon_set,
    left_shift_codon,
    # CyclesAnalysis.jl
    display_cycles,
    # GraphUtils.jl
    construct_graph!,
    create_vertices!,
    connect_edges!,
    add_vertex_by_label!,
    add_edge_by_label!,
    connect_edge_by_label!,
    create_shifted_graph,
    # PlotUtils.jl
    show_graph,
    show_multiple_graphs


"""
    demo_function(string::String)

A demo function that takes a string as input and turns it into all capitals.
# Example
```jldoctest
julia> demo_function("hello")
"HELLO"
```
"""

function demo_function(string::String)
    return uppercase(string)
end

end
