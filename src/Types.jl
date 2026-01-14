using BioSequences
using Graphs

# ---------------------------------------------- STRUCTS ----------------------------------------------


# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    graph::Graphs.SimpleDiGraph # directed graph
    codon_set::Vector{LongDNA{4}} # codon set represented in the graph
    all_vertex_labels::Vector{String} # all vertice labels
    base_vertex_labels::Vector{String} # base vertex labels
    added_vertex_labels::Vector{String} # added vertex labels
    all_edge_labels::Vector{Tuple{String, String}} # all edge labels
    base_edge_labels::Vector{Tuple{String, String}} # base edge labels
    added_edge_labels::Vector{Tuple{String, String}} # added edge labels
    vertex_index::Dict{String, Int} # from vertice label to vertice index ("AA" => 3)
    plot_title::String # title of the plot

    """
        CodonGraphData(codon_set::Vector{LongDNA{4}}; plot_title::String = "")

    Create a `CodonGraphData` instance with an empty graph and initialized fields.

    # Arguments
    - `codon_set::Vector{LongDNA{4}}`: Codon set for the graph.

    # Keyword Arguments
    - `plot_title::String`: Title for plotting (default: empty).

    # Returns
    - `CodonGraphData`: New graph data instance.

    # Throws
    - `ArgumentError`: If `codon_set` is empty.
    - `ArgumentError`: If any codon in `codon_set` does not have length 3.
    - `ArgumentError`: If `codon_set` contains duplicate codons.

    # Example
    ```julia
    codon_set = LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"])
    CodonGraphData(codon_set; plot_title = "Example")
    ```
    """
    # inner constructor to initialize empty graph, vectors and dicts
    function CodonGraphData(codon_set::Vector{LongDNA{4}}; plot_title::String = "")
        # do not allow empty codon sets
        isempty(codon_set) && throw(ArgumentError("Codon set cannot be empty!"))
        # each codon must have length 3
        any(length(codon) != 3 for codon in codon_set) &&
            throw(ArgumentError("All codons in codon set must have length 3!"))
        # do not allow duplicate codons in codon set
        length(codon_set) != length(Set(codon_set)) &&
            throw(ArgumentError("Codon set cannot contain duplicate codons!"))
        # only allow codons with valid DNA bases
        valid_bases = Set((DNA_A, DNA_C, DNA_G, DNA_T))
        any(any(base ∉ valid_bases for base in codon) for codon in codon_set) &&
            throw(ArgumentError("Codon set contains invalid DNA bases! Only A, C, G, T are allowed."))


        return new(
            Graphs.SimpleDiGraph(0),
            codon_set,
            String[],
            String[],
            String[],
            Tuple{String, String}[],
            Tuple{String, String}[],
            Tuple{String, String}[],
            Dict{String, Int}(),
            plot_title,
        )
    end
end
