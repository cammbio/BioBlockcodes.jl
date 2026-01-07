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
    is_codon_graphs_identical,
    _has_vertex_label,
    _has_edge_label,
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
    display_all_cycles,
    display_duplicate_cycles,
    get_cycle_count,
    get_cycles_of_length,
    get_cycle_count_of_length,
    get_cycles_from_vertex,
    get_cycles_difference,
    check_cycle_count_by_length,
    get_max_cycle_length,
    # GraphUtils.jl
    construct_graph_data!,
    _add_vertices_by_codon_set!,
    _add_edges_by_codon_set!,
    add_vertex_by_label!,
    add_edge_by_label!,
    _connect_edge_by_label!,
    create_shifted_graph,
    _check_data_fields_empty,
    # PlotUtils.jl
    show_graph,
    show_multiple_graphs,
    # WriterUtils.jl
    println_to_file
"""
    demo_function(string::String)

A demo function that takes a string as input and turns it into all capitals.

# Arguments
- `string::String`: Input string.

# Keyword Arguments
- None.

# Returns
- `String`: Uppercased string.

# Throws
- None.

# Example
```julia
demo_function("hello")
```
"""

function demo_function(string::String)
    return uppercase(string)
end

end
