using BioSequences
using BioBlockcodes
using Graphs
using Test


@testset "constructor guards" begin
    @testset "empty guard" begin
        codon_set = Vector{LongDNA{4}}()
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    @testset "length guard" begin
        codon_set = LongDNA{4}.(["CG", "GCA"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    @testset "duplicate guard" begin
        codon_set = LongDNA{4}.(["CGT", "GCA", "TTA", "CGT"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end


    @testset "valid bases guard" begin
        codon_set = LongDNA{4}.(["CGC", "GMA", "AGA"])
        @test_throws ArgumentError CodonGraphData(codon_set)
    end
end


@testset "field construction" begin
    @testset "graph field" begin
        codon_set = LongDNA{4}.(["CGT", "GCA"])
        cgd = CodonGraphData(codon_set)
        @test nv(cgd.graph) == 8
        @test ne(cgd.graph) == 4
    end


    @testset "codon_set field" begin
        codon_set = LongDNA{4}.(["GGT", "ACA"])
        cgd = CodonGraphData(codon_set)
        @test cgd.codon_set == codon_set
    end


    @testset "vert_labels field" begin
        codon_set = LongDNA{4}.(["ATG", "TTA"])
        cgd = CodonGraphData(codon_set)
        expected_vert_labels = ["A", "G", "T", "AT", "TA", "TG", "TT"]
        @test cgd.vert_labels == expected_vert_labels
    end


    @testset "edge_labels field" begin
        codon_set = LongDNA{4}.(["CCG", "GAC"])
        cgd = CodonGraphData(codon_set)
        expected_edge_labels = [("C", "CG"), ("CC", "G"), ("G", "AC"), ("GA", "C")]
        @test cgd.edge_labels == expected_edge_labels
    end


    @testset "vert_idxs field" begin
        codon_set = LongDNA{4}.(["TCG", "GGT"])
        cgd = CodonGraphData(codon_set)
        expected_vert_idxs = Dict("T" => 2, "GG" => 4, "CG" => 3, "G" => 1, "GT" => 5, "TC" => 6)
        @test cgd.vert_idxs == expected_vert_idxs
    end


    @testset "graph_title field" begin
        codon_set = LongDNA{4}.(["AAT", "AGC"])
        graph_title = "Test Graph"
        cgd = CodonGraphData(codon_set; graph_title = graph_title)
        @test cgd.graph_title == graph_title
    end
end