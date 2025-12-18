using BioSequences
using Graphs
using Test

println("Running Types tests...")


@testset "CodonGraphData constructor" begin
    # happy path
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAG", "CTT", "GGA"])
        plot_title = "My Codon Graph"
        data = CodonGraphData(codon_set; plot_title = plot_title)

        @test isa(data, CodonGraphData)
        @test data.codon_set == codon_set
        @test nv(data.graph) == 0
        @test ne(data.graph) == 0
        @test isempty(data.all_vertex_labels)
        @test isempty(data.base_vertex_labels)
        @test isempty(data.added_vertex_labels)
        @test isempty(data.all_edge_labels)
        @test isempty(data.base_edge_labels)
        @test isempty(data.added_edge_labels)
        @test isempty(data.vertex_index)
        @test data.plot_title == plot_title
    end

    # empty codon set (should throw ArgumentError)
    @testset "empty codon set" begin
        codon_set = LongDNA{4}[]
        @test_throws ArgumentError CodonGraphData(codon_set)
    end

    # codon set containing codon with length not equal to 3
    @testset "invalid codon length" begin
        codon_set = LongDNA{4}.(["AAT", "GT"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end

    # codon set containing duplicate codons
    @testset "duplicate codons" begin
        codon_set = LongDNA{4}.(["AAT", "GTC", "AAT"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------