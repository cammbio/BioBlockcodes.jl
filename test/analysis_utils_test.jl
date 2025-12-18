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

    # test graph with single vertex and no edges, returns ArgumentError
    @testset "single vertex no edges" begin
        graph = SimpleDiGraph(1)
        @test_throws ArgumentError is_circular(graph)
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

    # test graph where cycle does not start at vertex 1, returns false
    @testset "cycle not starting at vertex 1" begin
        graph = SimpleDiGraph(4)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 2)
        @test !is_circular(graph)
    end

    # test graph with disconnected acyclic components, returns true
    @testset "disconnected acyclic components" begin
        graph = SimpleDiGraph(6)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 6)
        @test is_circular(graph)
    end

    # test graph with multiple cycles, returns false
    @testset "multiple cycles" begin
        graph = SimpleDiGraph(4)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 1)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 3)
        add_edge!(graph, 1, 3)
        add_edge!(graph, 4, 2)
        @test !is_circular(graph)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "is_comma_free" begin
    # test graph with exactly one path, path shorter than 2, returns true
    @testset "one path shorter than 2" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 2)
        @test is_comma_free(graph)
    end

    # test graph with exactly one path, path equal to 2, returns true
    @testset "one path equal to 2" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        @test is_comma_free(graph)
    end

    # test graph with exactly one path, path longer than 2, returns false
    @testset "one path longer than 2" begin
        graph = SimpleDiGraph(4)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths and no disconnected components, all shorter than 2, returns true
    @testset "multiple paths and no disconnected components, all shorter than 2" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 3, 1)
        @test is_comma_free(graph)
    end

    # test graph with multiple paths and no disconnected components, all equal to 2, returns true
    @testset "multiple paths and no disconnected components, all equal to 2" begin
        graph = SimpleDiGraph(5)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 1, 4)
        add_edge!(graph, 4, 5)
        @test is_comma_free(graph)
    end

    # test graph with multiple paths and no disconnected components, all longer than 2, returns false
    @testset "multiple paths and no disconnected components, all longer than 2" begin
        graph = SimpleDiGraph(8)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 1, 5)
        add_edge!(graph, 5, 6)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths and no disconnected components, one shorter, one equal to, one longer than 2, returns false
    @testset "multiple paths and no disconnected components, one shorter, one equal to, one longer than 2" begin
        graph = SimpleDiGraph(8)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 1, 4)
        add_edge!(graph, 4, 3)
        add_edge!(graph, 1, 5)
        add_edge!(graph, 5, 6)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        add_edge!(graph, 8, 3)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths and disconnected components, all paths shorter than 2, returns true
    @testset "multiple paths and disconnected components, all shorter than 2" begin
        graph = SimpleDiGraph(6)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 5, 6)
        @test is_comma_free(graph)
    end

    # test graph with multiple paths and disconnected components, all paths equal to 2, returns true
    @testset "multiple paths and disconnected components, all equal to 2" begin
        graph = SimpleDiGraph(9)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 5)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 9)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        @test is_comma_free(graph)
    end

    # test graph with multiple paths and disconnected components, all paths longer than 2, returns false
    @testset "multiple paths and disconnected components, all longer than 2" begin
        graph = SimpleDiGraph(15)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        add_edge!(graph, 8, 9)
        add_edge!(graph, 10, 11)
        add_edge!(graph, 11, 12)
        add_edge!(graph, 12, 13)
        add_edge!(graph, 13, 14)
        add_edge!(graph, 14, 15)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths and disconnected components, one shorter, one equal to, one longer than 2, returns false
    @testset "multiple paths and disconnected components, one shorter, one equal to, one longer than 2" begin
        graph = SimpleDiGraph(12)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        add_edge!(graph, 8, 9)
        add_edge!(graph, 9, 10)
        add_edge!(graph, 10, 11)
        add_edge!(graph, 11, 12)
        @test !is_comma_free(graph)
    end

    # test graph with one path which is a self-loop, returns false
    @testset "one path, self-loop" begin
        graph = SimpleDiGraph(1)
        add_edge!(graph, 1, 1)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, one is a self-loop, returns false
    @testset "multiple paths, one self-loop" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 1)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, some are self-loops, returns false
    @testset "multiple paths, some self-loops" begin
        graph = SimpleDiGraph(3)
        add_edge!(graph, 1, 1)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 3)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, all are self-loops, returns false
    @testset "multiple paths, all self-loops" begin
        graph = SimpleDiGraph(4)
        add_edge!(graph, 1, 1)
        add_edge!(graph, 2, 2)
        add_edge!(graph, 3, 3)
        add_edge!(graph, 4, 4)
        @test !is_comma_free(graph)
    end

    # test graph with one path which is a normal loop between two vertices, returns false
    @testset "one path, normal loop between two vertices" begin
        graph = SimpleDiGraph(2)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 1)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, one is a normal loop between two vertices, returns false
    @testset "multiple paths, one normal loop between two vertices" begin
        graph = SimpleDiGraph(5)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 1)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 4, 5)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, some are normal loops between two vertices, returns false
    @testset "multiple paths, some normal loops between two vertices" begin
        graph = SimpleDiGraph(6)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 1)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 6)
        add_edge!(graph, 6, 5)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, all are normal loops between two vertices, returns false
    @testset "multiple paths, all normal loops between two vertices" begin
        graph = SimpleDiGraph(4)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 1)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 2)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 3)
        @test !is_comma_free(graph)
    end

    # test graph with one path which is a normal loop between five vertices, returns false
    @testset "one path, normal loop between five vertices" begin
        graph = SimpleDiGraph(5)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 1)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, one is a normal loop between five vertices, returns false
    @testset "multiple paths, one normal loop between five vertices" begin
        graph = SimpleDiGraph(8)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 1)
        add_edge!(graph, 1, 6)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, some are normal loops between five vertices, returns false
    @testset "multiple paths, some normal loops between three to five vertices" begin
        graph = SimpleDiGraph(12)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 1)
        add_edge!(graph, 2, 6)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        add_edge!(graph, 8, 2)
        add_edge!(graph, 2, 6)
        add_edge!(graph, 9, 10)
        add_edge!(graph, 10, 6)
        add_edge!(graph, 6, 9)
        add_edge!(graph, 1, 11)
        add_edge!(graph, 11, 12)
        @test !is_comma_free(graph)
    end

    # test graph with multiple paths, all are normal loops between five vertices, returns false
    @testset "multiple paths, all normal loops between three to five vertices" begin
        graph = SimpleDiGraph(8)
        add_edge!(graph, 1, 2)
        add_edge!(graph, 2, 3)
        add_edge!(graph, 3, 4)
        add_edge!(graph, 4, 5)
        add_edge!(graph, 5, 1)
        add_edge!(graph, 1, 6)
        add_edge!(graph, 6, 7)
        add_edge!(graph, 7, 8)
        add_edge!(graph, 8, 1)
        add_edge!(graph, 2, 8)
        @test !is_comma_free(graph)
    end

    # test graph with no vertices, throws ArgumentError
    @testset "empty graph" begin
        graph = SimpleDiGraph(0)
        @test_throws ArgumentError is_comma_free(graph)
    end

    # test graph with vertices but no edges, throws ArgumentError
    @testset "graph with no edges" begin
        graph = SimpleDiGraph(3)
        @test_throws ArgumentError is_comma_free(graph)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "is_self_complementary" begin
    # test self-complementary codon set, returns true
    @testset "self-complementary codon set" begin
        codon_set = LongDNA{4}.(["ATG", "CAT", "CTT", "AAG"])
        data = CodonGraphData(codon_set)
        construct_graph!(data; show_debug = false, show_plot = false)
        @test is_self_complementary(data; show_debug = false, show_plot = false)
    end

    # test non-self-complementary codon set, returns false
    @testset "non-self-complementary codon set" begin
        codon_set = LongDNA{4}.(["AGT", "TGA", "CAA", "TGT", "GGA"])
        data = CodonGraphData(codon_set)
        construct_graph!(data; show_debug = false, show_plot = false)
        @test !is_self_complementary(data; show_debug = false, show_plot = false)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "is_graphs_identical" begin

end

# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
