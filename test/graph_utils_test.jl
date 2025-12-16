using BioSequences
using Graphs
using Test

println("Running GraphUtils tests...")

@testset "construct_graph!" begin
    @testset "properties" begin
        codon_set = LongDNA{4}.(["GTA", "GTT", "GCA", "GCG"])
        data = CodonGraphData(codon_set)
        construct_graph!(data; show_plot = false, show_debug = false)

        # test properties of the constructed graph and data
        @test nv(data.graph) == 9
        @test ne(data.graph) == 8
        @test data.codon_set == codon_set
        @test data.all_vertex_labels == ["A", "G", "T", "CA", "CG", "GC", "GT", "TA", "TT"]
        @test data.base_vertex_labels == ["A", "G", "T", "CA", "CG", "GC", "GT", "TA", "TT"]
        @test data.added_vertex_labels == String[]
        @test data.all_edge_labels == [
            ("G", "TA"),
            ("GT", "A"),
            ("G", "TT"),
            ("GT", "T"),
            ("G", "CA"),
            ("GC", "A"),
            ("G", "CG"),
            ("GC", "G"),
        ]
        @test data.base_edge_labels == [
            ("G", "TA"),
            ("GT", "A"),
            ("G", "TT"),
            ("GT", "T"),
            ("G", "CA"),
            ("GC", "A"),
            ("G", "CG"),
            ("GC", "G"),
        ]
        @test data.added_edge_labels == Tuple{String, String}[]
        @test data.vertex_index == Dict(
            "A" => 1,
            "G" => 2,
            "T" => 3,
            "CA" => 4,
            "CG" => 5,
            "GC" => 6,
            "GT" => 7,
            "TA" => 8,
            "TT" => 9,
        )

        @test data.base_vertex_labels[data.vertex_index["GC"]] == "GC"
        @test data.base_vertex_labels[data.vertex_index["TT"]] != "A"
    end
end