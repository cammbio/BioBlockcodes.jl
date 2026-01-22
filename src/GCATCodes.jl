module GCATCodes

# load other source files
include("WriteUtils.jl")
include("Types.jl")
include("CyclesAnalysis.jl") # needs Types.jl
include("CodonUtils.jl") # needs Types.jl
include("PlotUtils.jl") # needs Types.jl
include("CodonGraphAnalysis.jl") # needs Types.jl and CodonUtils.jl
include("GraphUtils.jl") # needs Types.jl and PlotUtils.jl
include("StrongC3Analysis.jl") # needs Types.jl, CodonUtils.jl and CodonGraphAnalysis.jl


# export types
export
# Types.jl
    CodonGraphData
#

# export functions
export
    # CodonGraphAnalysis.jl
    is_circular,
    is_c3,
    is_codon_graphs_identical,
    is_comma_free,
    is_self_complementary,
    is_strong_c3,
    # CodonUtils.jl
    get_codon_combinations_per_size,
    get_complemented_reversed_codon_set,
    get_complemented_codon_set,
    get_reversed_codon_set,
    get_complemented_codon,
    get_reversed_codon,
    get_complemented_base,
    left_shift_codon_set,
    left_shift_codon,
    # CyclesAnalysis.jl
    get_cycle_count,
    get_cycle_count_by_length,
    get_cycles_all,
    get_cycles_by_length,
    get_cycles_by_vertex_label,
    get_cycles_difference,
    get_duplicate_cycles,
    get_max_cycle_length,
    # GraphUtils.jl
    construct_graph_data!,
    add_vertex_by_label!,
    add_edge_by_label!,
    # PlotUtils.jl
    show_codon_graph,
    show_multiple_codon_graphs,
    # WriterUtils.jl
    line_to_codon_set,
    print_to_file,
    # strong C3 analysis
    process_strong_c3_combinations,
    process_strong_c3_combinations_by_combination_size,


    # temporary helper
    _expand_graph,
    _add_n2_n3n1_by_codon,
    _has_cycle_longer_than,
    _increment_codon_set_combination!,
    _save_strong_c3_checkpoint!,
    _load_strong_c3_checkpoint,
    _get_last_combination_indices_from_file,
    _get_next_combination,
    _get_last_line,
    _remove_empty_last_lines,
    _contains_codon_rotation
end