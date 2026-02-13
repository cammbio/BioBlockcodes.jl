using Test
using BioSequences
using GCATCodes
using Graphs


println("Running tests for Types.jl...")


# test constructor guards
@testset "constructor guards" begin
    # codon set cannot be empty
    @testset "empty guard" begin
        codon_set = Vector{LongDNA{4}}()
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    # codons in set must have length 3
    @testset "length guard" begin
        codon_set = LongDNA{4}.(["CG", "GCA"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    # codon set cannot contain duplicate codons
    @testset "duplicate guard" begin
        codon_set = LongDNA{4}.(["CGT", "GCA", "TTA", "CGT"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    # codon set cannot contain invalid DNA bases (only A, C, G, T allowed)
    @testset "valid bases guard" begin
        codon_set = LongDNA{4}.(["CGC", "GMA", "AGA"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end
end


# test that fields are constructed correctly
@testset "field construction" begin
    # graph has correct number of vertices and edges based on codon set
    @testset "graph field" begin
        codon_set = LongDNA{4}.(["CGT", "GCA"])
        data = CodonGraphData(codon_set)
        @test nv(data.graph) == 8
        @test ne(data.graph) == 4
    end


    # codon_set field is correctly set to input codon set
    @testset "codon_set field" begin
        codon_set = LongDNA{4}.(["GGT", "ACA"])
        data = CodonGraphData(codon_set)
        @test data.codon_set == codon_set
    end


    # vert_labels field is correctly constructed based on codon set
    @testset "vert_labels field" begin
        codon_set = LongDNA{4}.(["ATG", "TTA"])
        data = CodonGraphData(codon_set)
        expected_vert_labels = ["A", "G", "T", "AT", "TA", "TG", "TT"]
        @test data.vert_labels == expected_vert_labels
    end


    # edge_labels field is correctly constructed based on codon set
    @testset "edge_labels field" begin
        codon_set = LongDNA{4}.(["CCG", "GAC"])
        data = CodonGraphData(codon_set)
        expected_edge_labels = [("C", "CG"), ("CC", "G"), ("G", "AC"), ("GA", "C")]
        @test data.edge_labels == expected_edge_labels
    end


    # vert_idxs field is correctly constructed as a mapping from vert_labels to vertex IDs
    @testset "vert_idxs field" begin
        codon_set = LongDNA{4}.(["TCG", "GGT"])
        data = CodonGraphData(codon_set)
        expected_vert_idxs = Dict("T" => 2, "GG" => 4, "CG" => 3, "G" => 1, "GT" => 5, "TC" => 6)
        @test data.vert_idxs == expected_vert_idxs
    end


    # graph_title field is correctly set to input graph_title
    @testset "graph_title field" begin
        codon_set = LongDNA{4}.(["AAT", "AGC"])
        graph_title = "Test Graph"
        data = CodonGraphData(codon_set; graph_title = graph_title)
        @test data.graph_title == graph_title
    end
end