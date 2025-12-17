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

    # test empty graph, throws ArgumentError
    @testset "empty graph" begin
        graph = SimpleDiGraph(0)
        @test_throws ArgumentError is_circular(graph)
    end

    # test graph with self-loop, returns false
    @testset "self-loop graph" begin
        graph = SimpleDiGraph(2)
        add_edge!(graph, 1, 1)
        add_edge!(graph, 1, 2)
        @test !is_circular(graph)
    end

    # test graph with single vertex and no edges, returns true
    @testset "single vertex no edges" begin
        graph = SimpleDiGraph(1)
        @test is_circular(graph)
    end

    # test graph with disconnected components, one cyclic and one acyclic, returns false
    @testset "disconnected components" begin
        graph = SimpleDiGraph(5)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 4)
        @test !is_circular(graph)
    end
end

