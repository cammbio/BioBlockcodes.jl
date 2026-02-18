using Test
using BioSequences
using GCATCodes
using Graphs


println("Running tests for CodonGraphAnalysis.jl...")


@testset "is_c3" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        data = CodonGraphData(codon_set)
        @test is_c3(data) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "GAA"])
        data = CodonGraphData(codon_set)
        @test is_c3(data) == false
    end


    @testset "alpha_1 not circular" begin
        codon_set = LongDNA{4}.(["AAC", "CCA"])
        data = CodonGraphData(codon_set)
        @test is_c3(data) == false
    end


    @testset "alpha_2 not circular" begin
        codon_set = LongDNA{4}.(["TGG", "GTT"])
        data = CodonGraphData(codon_set)
        @test is_c3(data) == false
    end
end


@testset "is_circular" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        data = CodonGraphData(codon_set)
        @test is_circular(data) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "GAA"])
        data = CodonGraphData(codon_set)
        @test is_circular(data) == false
    end
end


@testset "is_comma_free" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        data = CodonGraphData(codon_set)
        @test is_comma_free(data) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "CTG", "TGA", "TTC"])
        data = CodonGraphData(codon_set)
        @test is_comma_free(data) == false
    end
end


@testset "is_self_complementary" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "GTT"])
        data = CodonGraphData(codon_set)
        @test is_self_complementary(data) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG"])
        data = CodonGraphData(codon_set)
        @test is_self_complementary(data) == false
    end
end


@testset "is_strong_c3" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        data = CodonGraphData(codon_set)
        @test is_strong_c3(data) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "CTG", "TGA", "TTC"])
        data = CodonGraphData(codon_set)
        @test is_strong_c3(data) == false
    end


    @testset "is_c3 true but not is_strong_c3" begin
        codon_set = LongDNA{4}.(["AAC", "ACC"])
        data = CodonGraphData(codon_set)
        @test is_strong_c3(data) == false
    end
end


@testset "_add_edge_by_label!" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["AGT", "GCC"])


    @testset "happy path" begin
        data = CodonGraphData(codon_set)
        edge_label = ("T", "CC")
        src_label, dst_label = edge_label
        @test GCATCodes._add_edge_by_label!(data, edge_label) == true
        @test (edge_label in data.edge_labels) == true
        @test (has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])) == true
        @test ne(data.graph) == 5
    end


    @testset "sad path" begin
        data = CodonGraphData(codon_set)
        edge_label = ("A", "GT")
        src_label, dst_label = edge_label
        @test GCATCodes._add_edge_by_label!(data, edge_label) == false
        @test (edge_label in data.edge_labels) == true
        @test (has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])) == true
        @test ne(data.graph) == 4
    end


    @testset "src_label not in vert_labels" begin
        data = CodonGraphData(codon_set)
        edge_label = ("AT", "C")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(data, edge_label)
        @test (edge_label in data.edge_labels) == false
        @test_throws KeyError has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])
    end


    @testset "dst_label not in vert_labels" begin
        data = CodonGraphData(codon_set)
        edge_label = ("A", "TT")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(data, edge_label)
        @test (edge_label in data.edge_labels) == false
        @test_throws KeyError has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])
    end


    @testset "src_label not in vert_idxs" begin
        data = CodonGraphData(codon_set)
        edge_label = ("CG", "T")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(data, edge_label)
        @test (edge_label in data.edge_labels) == false
        @test_throws KeyError has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])
    end


    @testset "dst_label not in vert_idxs" begin
        data = CodonGraphData(codon_set)
        edge_label = ("G", "CT")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(data, edge_label)
        @test (edge_label in data.edge_labels) == false
        @test_throws KeyError has_edge(data.graph, data.vert_idxs[src_label], data.vert_idxs[dst_label])
    end
end


@testset "_add_vert_by_label!" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CCG", "GTA"])


    @testset "happy path" begin
        data = CodonGraphData(codon_set)
        vert_label = "AC"
        @test GCATCodes._add_vert_by_label!(data, vert_label) == true
        @test (vert_label in data.vert_labels) == true
        @test haskey(data.vert_idxs, vert_label) == true
        @test nv(data.graph) == 8
    end


    @testset "sad path" begin
        data = CodonGraphData(codon_set)
        vert_label = "CC"
        @test GCATCodes._add_vert_by_label!(data, vert_label) == false
        @test (vert_label in data.vert_labels) == true
        @test haskey(data.vert_idxs, vert_label) == true
        @test nv(data.graph) == 7
    end
end


@testset "_dfs_depth" begin
    depth = 0
    max_depth = 4
    start_vert = 1


    @testset "normal path length < max_depth" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == false
    end


    @testset "normal path length == max_depth" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == false
    end


    @testset "normal path length > max_depth" begin
        g = SimpleDiGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == true
    end


    @testset "loop path length < max_depth" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        # @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == false
    end


    @testset "loop path length == max_depth" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 1)
        # @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == false
    end


    @testset "loop path length > max_depth" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 1)
        @test GCATCodes._dfs_depth(g, start_vert, depth, max_depth) == true
    end
end


@testset "_expand_graph!" begin
    @testset "validate edges" begin
        codon_set = LongDNA{4}.(["AGA", "GAC", "TGG"])
        data = CodonGraphData(codon_set)
        _expand_graph!(data)
        @test ne(data.graph) == 12
        @test data.edge_labels == [
            ("A", "GA"),
            ("AG", "A"),
            ("G", "AC"),
            ("GA", "C"),
            ("T", "GG"),
            ("TG", "G"),
            ("G", "AA"),
            ("AA", "G"),
            ("A", "CG"),
            ("CG", "A"),
            ("G", "GT"),
            ("GT", "G"),
        ]
    end


    @testset "validate vertices" begin
        codon_set = LongDNA{4}.(["AAC", "ATG", "TGC"])
        data = CodonGraphData(codon_set)
        _expand_graph!(data)
        @test nv(data.graph) == 12
        @test data.vert_labels == ["A", "C", "G", "T", "AA", "AC", "AT", "GC", "TG", "CA", "GA", "CT"]
    end
end


@testset "_has_cycle_longer_than" begin
    @testset "no cycles" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        @test GCATCodes._has_cycle_longer_than(g, 0) == false
        @test GCATCodes._has_cycle_longer_than(g, 1) == false
        @test GCATCodes._has_cycle_longer_than(g, 4) == false
        @test GCATCodes._has_cycle_longer_than(g, 7) == false
    end


    @testset "self-loop" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 1)
        add_edge!(g, 2, 2)
        @test GCATCodes._has_cycle_longer_than(g, 0) == true
        @test GCATCodes._has_cycle_longer_than(g, 1) == false
        @test GCATCodes._has_cycle_longer_than(g, 4) == false
        @test GCATCodes._has_cycle_longer_than(g, 7) == false
    end


    @testset "normal loop of length 2" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        @test GCATCodes._has_cycle_longer_than(g, 0) == true
        @test GCATCodes._has_cycle_longer_than(g, 1) == true
        @test GCATCodes._has_cycle_longer_than(g, 2) == false
        @test GCATCodes._has_cycle_longer_than(g, 5) == false
    end


    @testset "normal loop of length 5" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 1)
        @test GCATCodes._has_cycle_longer_than(g, 0) == true
        @test GCATCodes._has_cycle_longer_than(g, 1) == true
        @test GCATCodes._has_cycle_longer_than(g, 4) == true
        @test GCATCodes._has_cycle_longer_than(g, 7) == false
    end


    @testset "a self-loop and a normal loop" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 1)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        @test GCATCodes._has_cycle_longer_than(g, 0) == true
        @test GCATCodes._has_cycle_longer_than(g, 1) == true
        @test GCATCodes._has_cycle_longer_than(g, 4) == false
        @test GCATCodes._has_cycle_longer_than(g, 7) == false
    end


    @testset "invalid no vertices and edges" begin
        g = SimpleDiGraph(0)
        @test_throws ArgumentError _has_cycle_longer_than(g, 1)
    end


    @testset "invalid vertices but no edges" begin
        g = SimpleDiGraph(3)
        @test_throws ArgumentError _has_cycle_longer_than(g, 1)
    end


    @testset "invalid max_length argument" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        @test_throws ArgumentError _has_cycle_longer_than(g, -2)
    end
end


@testset "_is_codon_graphs_equal" begin
    @testset "happy path" begin
        codon_set_1 = LongDNA{4}.(["AAC", "ATG"])
        codon_set_2 = LongDNA{4}.(["AAC", "ATG"])
        data_1 = CodonGraphData(codon_set_1)
        data_2 = CodonGraphData(codon_set_2)
        @test GCATCodes._is_codon_graphs_equal(data_1, data_2) == true
    end


    @testset "sad path" begin
        codon_set_1 = LongDNA{4}.(["AGT", "GCA"])
        codon_set_2 = LongDNA{4}.(["AGT", "CGA"])
        data_1 = CodonGraphData(codon_set_1)
        data_2 = CodonGraphData(codon_set_2)
        @test GCATCodes._is_codon_graphs_equal(data_1, data_2) == false
    end


    @testset "same vertices but different edges" begin
        codon_set = LongDNA{4}.(["GAA", "GCA"])
        data_1 = CodonGraphData(codon_set)
        data_2 = CodonGraphData(codon_set)
        edge_label = ("A", "CA")
        GCATCodes._add_edge_by_label!(data_1, edge_label)
        @test GCATCodes._is_codon_graphs_equal(data_1, data_2) == false
    end


    @testset "same edges but different vertices" begin
        codon_set_1 = LongDNA{4}.(["AGT", "GCA"])
        codon_set_2 = LongDNA{4}.(["AGT", "GCA"])
        data_1 = CodonGraphData(codon_set_1)
        data_2 = CodonGraphData(codon_set_2)
        GCATCodes._add_vert_by_label!(data_1, "CT")
        @test GCATCodes._is_codon_graphs_equal(data_1, data_2) == false
    end
end


@testset "_has_edge_label" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CGT", "TTA"])
    data = CodonGraphData(codon_set)

    @testset "happy path" begin
        @test GCATCodes._has_edge_label(data, "C", "GT") == true
    end


    @testset "sad path" begin
        @test GCATCodes._has_edge_label(data, "C", "TA") == false
    end


    @testset "src_label not in vert_labels" begin
        @test_throws ArgumentError GCATCodes._has_edge_label(data, "GA", "T")
    end


    @testset "dst_label not in vert_labels" begin
        @test_throws ArgumentError GCATCodes._has_edge_label(data, "T", "AT")
    end


    @testset "src_label not in vert_idxs" begin
        @test_throws ArgumentError GCATCodes._has_edge_label(data, "G", "AA")
    end


    @testset "dst_label not in vert_idxs" begin
        @test_throws ArgumentError GCATCodes._has_edge_label(data, "TT", "G")
    end
end


@testset "_has_vert_label" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CTT", "ACT"])
    data = CodonGraphData(codon_set)


    @testset "happy path" begin
        @test GCATCodes._has_vert_label(data, "CT") == true
    end


    @testset "sad path" begin
        @test GCATCodes._has_vert_label(data, "AG") == false
    end
end