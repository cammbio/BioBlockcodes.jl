module BioBlockcodes
# ----------------------------------------------- PACKAGES ----------------------------------------------
using Base.Threads: @spawn
using BioSequences
using CairoMakie
using GraphMakie
using Graphs
using NetworkLayout
# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------
const ALL_CODONS =
    LongDNA{
        4,
    }.([
        "AAC",
        "AAG",
        "AAT",
        "ACA",
        "ACC",
        "ACG",
        "ACT",
        "AGA",
        "AGC",
        "AGG",
        "AGT",
        "ATA",
        "ATC",
        "ATG",
        "ATT",
        "CAA",
        "CAC",
        "CAG",
        "CAT",
        "CCA",
        "CCG",
        "CCT",
        "CGA",
        "CGC",
        "CGG",
        "CGT",
        "CTA",
        "CTC",
        "CTG",
        "CTT",
        "GAA",
        "GAC",
        "GAG",
        "GAT",
        "GCA",
        "GCC",
        "GCG",
        "GCT",
        "GGA",
        "GGC",
        "GGT",
        "GTA",
        "GTC",
        "GTG",
        "GTT",
        "TAA",
        "TAC",
        "TAG",
        "TAT",
        "TCA",
        "TCC",
        "TCG",
        "TCT",
        "TGA",
        "TGC",
        "TGG",
        "TGT",
        "TTA",
        "TTC",
        "TTG",
    ])

const ALLOWED_BASES_DNA = Set((DNA_A, DNA_C, DNA_G, DNA_T))
const ALLOWED_BASES_STR = Set(('A', 'C', 'G', 'T'))
const CODON_INDEX = Dict{LongDNA{4}, Int}(codon => i for (i, codon) in enumerate(ALL_CODONS))
# ---------------------------------------------- INCLUDES -----------------------------------------------
include("Types.jl")
include("Validation.jl")
include("InputOutputUtilities.jl")
include("CodonUtilities.jl") # needs Types.jl
include("CodonGraphAnalysis.jl") # needs Types.jl, CodonUtilities.jl and CodonGraphBuilder.jl
include("CycleAnalysis.jl") # needs Types.jl
include("CodonGraphPlotting.jl") # needs Types.jl
include("StrongC3.jl") # needs Types.jl, CodonUtilities.jl and CodonGraphAnalysis.jl
# ---------------------------------------------- EXPORTS ------------------------------------------------
# Types
export
# Types.jl
    CodonGraphData

# Functions
export
    # CodonGraphAnalysis.jl
    is_c3,
    is_circular,
    is_comma_free,
    is_self_complementary,
    is_strong_c3,
    # CodonGraphPlotting.jl
    plot_codon_graph,
    plot_multiple_codon_graphs,
    # CodonUtilites.jl
    get_comp_rev_codon_set,
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
    # InputOutputUtils.jl
    get_codon_set_from_line,
    codon_set_to_str,
    # StrongC3.jl
    calc_strong_c3_comb_by_size,

    _expand_graph!
end
