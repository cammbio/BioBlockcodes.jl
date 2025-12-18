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
        construct_graph!(data)
        @test is_self_complementary(data)
    end

    # test non-self-complementary codon set, returns false
    @testset "non-self-complementary codon set" begin
        codon_set = LongDNA{4}.(["AGT", "TGA", "CAA", "TGT", "GGA"])
        data = CodonGraphData(codon_set)
        construct_graph!(data)
        @test !is_self_complementary(data)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "is_graphs_identical" begin
    # test identical graphs from CodonGraphData objects, returns true
    @testset "identical graphs" begin
        data_1 = CodonGraphData(LongDNA{4}.(["ATA", "AAG", "CAC"]))
        construct_graph!(data_1)
        data_2 = CodonGraphData(LongDNA{4}.(["ATA", "AAG", "CAC"]))
        construct_graph!(data_2)
        @test is_graphs_identical(data_1, data_2)
    end

    # test identical graphs from CodonGraphData objects with same codon set but different order, returns true
    @testset "identical graphs with different codon order" begin
        data_1 = CodonGraphData(LongDNA{4}.(["GTC", "TTA", "AGC"]))
        construct_graph!(data_1)
        data_2 = CodonGraphData(LongDNA{4}.(["AGC", "GTC", "TTA"]))
        construct_graph!(data_2)
        @test is_graphs_identical(data_1, data_2)
    end

    # test non-identical graphs from CodonGraphData objects, returns false
    @testset "non-identical graphs" begin
        data_1 = CodonGraphData(LongDNA{4}.(["AGT", "ACG", "ATC"]))
        construct_graph!(data_1)
        data_2 = CodonGraphData(LongDNA{4}.(["AGT", "ACG", "ATT"]))
        construct_graph!(data_2)
        @test !is_graphs_identical(data_1, data_2)
    end

    # test method call symmetry, returns same result both ways
    @testset "method call symmetry" begin
        data_1 = CodonGraphData(LongDNA{4}.(["TAC", "GGA", "CTT"]))
        construct_graph!(data_1)
        data_2 = CodonGraphData(LongDNA{4}.(["TAC", "GGA", "CTT"]))
        construct_graph!(data_2)
        @test is_graphs_identical(data_1, data_2) == is_graphs_identical(data_2, data_1)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "has_vertex_label" begin
    # test existing vertex label, returns true
    @testset "existing vertex label" begin
        data = CodonGraphData(LongDNA{4}.(["ATG", "TAC", "GGA"]))
        construct_graph!(data)
        @test has_vertex_label(data.vertex_index, "AT")
    end

    # test non-existing vertex label, returns false
    @testset "non-existing vertex label" begin
        data = CodonGraphData(LongDNA{4}.(["CCT", "GGA", "TTA"]))
        construct_graph!(data)
        @test !has_vertex_label(data.vertex_index, "AA")
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "has_edge_label" begin
    # test existing edge label, returns true
    @testset "existing edge label" begin
        data = CodonGraphData(LongDNA{4}.(["ATG", "TAC", "GGA"]))
        construct_graph!(data)
        @test has_edge_label(data, "AT", "G")
    end

    # test non-existing edge label, returns false
    @testset "non-existing edge label" begin
        data = CodonGraphData(LongDNA{4}.(["CCT", "GGA", "TTA"]))
        construct_graph!(data)
        @test !has_edge_label(data, "GG", "T")
    end

    # test edge label with non-existing vertices, returns false
    @testset "edge label with non-existing vertices" begin
        data = CodonGraphData(LongDNA{4}.(["CCT", "GGA", "TTA"]))
        construct_graph!(data)
        @test !has_edge_label(data, "CG", "CT")
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
@testset "is_c3" begin
    # test codon set that is C3, returns true
    @testset "C3 codon set" begin
        codon_set = LongDNA{4}.(["ATT", "TAC", "GGA", "CCT"])
        data = CodonGraphData(codon_set)
        construct_graph!(data)
        @test is_c3(data)
    end

    # test codon set that is not circular aka. not C3, returns false
    @testset "non-C3 codon set" begin
        codon_set = LongDNA{4}.(["ATG", "TAC", "GGA", "CCT", "AAA"])
        data = CodonGraphData(codon_set)
        construct_graph!(data)
        @test !is_c3(data)
    end

    # test codon set that is circular but not C3 aka. shifted graph is not circular, returns false
    @testset "circular but non-C3 codon set" begin
        codon_set = LongDNA{4}.(["AGT", "TAT", "CCT", "GAG", "AAC", "AAT", "GAT", "CCA"])
        data = CodonGraphData(codon_set)
        construct_graph!(data)
        @test !is_c3(data)
    end
end


# ------------------------------------ NEXT FUNCTION -------------------------------------------------------
