module GCATCodes

# load other source files
include("Types.jl")
include("CodonUtils.jl") # needs Types.jl
include("AnalysisUtils.jl") # needs Types.jl
include("PlotUtils.jl") # needs Types.jl
include("GraphUtils.jl") # needs Types.jl and PlotUtils.jl

# export relevant types
export
# Types.jl
    CodonGraphData
#

# export relevant functions
export
    # AnalysisUtils.jl
    add_vertex_by_label!,
    add_edge_by_label!,
    connect_edge_by_label!,
    display_cycles,
    has_edge_label,
    is_c3,
    is_circular,
    is_comma_free,
    is_graphs_identical,
    is_self_complementary,
    left_shift_codon_set,
    left_shift_codon,
    # CodonUtils.jl
    get_complement_reversed_codon_set,
    get_reverse_codon_set,
    # GCATCodes.jl
    # GraphUtils.jl
    construct_graph!,
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
