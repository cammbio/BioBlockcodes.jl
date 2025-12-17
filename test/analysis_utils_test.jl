using BioSequences
using Graphs
using Test

println("Running AnalysisUtils tests...")

@testset "is_circular" begin
    # test acyclic graph, returns true
    @testset "acyclic graph" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        @test is_circular(graph)
    end

    # test cyclic graph, returns false
    @testset "cyclic graph" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 1)
        @test !is_circular(graph)
    end
end

