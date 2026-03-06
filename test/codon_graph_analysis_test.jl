using BioSequences
using GCATCodes
using Graphs
using Test


@testset "is_c3" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        cgd = CodonGraphData(codon_set)
        @test is_c3(cgd) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "GAA"])
        cgd = CodonGraphData(codon_set)
        @test is_c3(cgd) == false
    end


    @testset "alpha_1 not circular" begin
        codon_set = LongDNA{4}.(["AAC", "CCA"])
        cgd = CodonGraphData(codon_set)
        @test is_c3(cgd) == false
    end


    @testset "alpha_2 not circular" begin
        codon_set = LongDNA{4}.(["TGG", "GTT"])
        cgd = CodonGraphData(codon_set)
        @test is_c3(cgd) == false
    end
end


@testset "is_circular" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        cgd = CodonGraphData(codon_set)
        @test is_circular(cgd) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "GAA"])
        cgd = CodonGraphData(codon_set)
        @test is_circular(cgd) == false
    end
end


@testset "is_comma_free" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        cgd = CodonGraphData(codon_set)
        @test is_comma_free(cgd) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "CTG", "TGA", "TTC"])
        cgd = CodonGraphData(codon_set)
        @test is_comma_free(cgd) == false
    end
end


@testset "is_self_complementary" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "GTT"])
        cgd = CodonGraphData(codon_set)
        @test is_self_complementary(cgd) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG"])
        cgd = CodonGraphData(codon_set)
        @test is_self_complementary(cgd) == false
    end
end


@testset "is_strong_c3" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAC", "ATG"])
        cgd = CodonGraphData(codon_set)
        @test is_strong_c3(cgd) == true
    end


    @testset "sad path" begin
        codon_set = LongDNA{4}.(["AGA", "AAG", "CTG", "TGA", "TTC"])
        cgd = CodonGraphData(codon_set)
        @test is_strong_c3(cgd) == false
    end


    @testset "is_c3 true but not is_strong_c3" begin
        codon_set = LongDNA{4}.(["AAC", "ACC"])
        cgd = CodonGraphData(codon_set)
        @test is_strong_c3(cgd) == false
    end
end


@testset "_add_edge_by_label!" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["AGT", "GCC"])


    @testset "happy path" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("T", "CC")
        src_label, dst_label = edge_label
        @test GCATCodes._add_edge_by_label!(cgd, edge_label) == true
        @test (edge_label in cgd.edge_labels) == true
        @test (has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])) == true
        @test ne(cgd.graph) == 5
    end


    @testset "sad path" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("A", "GT")
        src_label, dst_label = edge_label
        @test GCATCodes._add_edge_by_label!(cgd, edge_label) == false
        @test (edge_label in cgd.edge_labels) == true
        @test (has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])) == true
        @test ne(cgd.graph) == 4
    end


    @testset "src_label not in vert_labels" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("AT", "C")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(cgd, edge_label)
        @test (edge_label in cgd.edge_labels) == false
        @test_throws KeyError has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])
    end


    @testset "dst_label not in vert_labels" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("A", "TT")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(cgd, edge_label)
        @test (edge_label in cgd.edge_labels) == false
        @test_throws KeyError has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])
    end


    @testset "src_label not in vert_idxs" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("CG", "T")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(cgd, edge_label)
        @test (edge_label in cgd.edge_labels) == false
        @test_throws KeyError has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])
    end


    @testset "dst_label not in vert_idxs" begin
        cgd = CodonGraphData(codon_set)
        edge_label = ("G", "CT")
        src_label, dst_label = edge_label
        @test_throws ArgumentError GCATCodes._add_edge_by_label!(cgd, edge_label)
        @test (edge_label in cgd.edge_labels) == false
        @test_throws KeyError has_edge(cgd.graph, cgd.vert_idxs[src_label], cgd.vert_idxs[dst_label])
    end
end


@testset "_add_vert_by_label!" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CCG", "GTA"])


    @testset "happy path" begin
        cgd = CodonGraphData(codon_set)
        vert_label = "AC"
        @test GCATCodes._add_vert_by_label!(cgd, vert_label) == true
        @test (vert_label in cgd.vert_labels) == true
        @test haskey(cgd.vert_idxs, vert_label) == true
        @test nv(cgd.graph) == 8
    end


    @testset "sad path" begin
        cgd = CodonGraphData(codon_set)
        vert_label = "CC"
        @test GCATCodes._add_vert_by_label!(cgd, vert_label) == false
        @test (vert_label in cgd.vert_labels) == true
        @test haskey(cgd.vert_idxs, vert_label) == true
        @test nv(cgd.graph) == 7
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
        cgd = CodonGraphData(codon_set)
        GCATCodes._expand_graph!(cgd)
        @test ne(cgd.graph) == 12
        @test cgd.edge_labels == [
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
        cgd = CodonGraphData(codon_set)
        GCATCodes._expand_graph!(cgd)
        @test nv(cgd.graph) == 12
        @test cgd.vert_labels == ["A", "C", "G", "T", "AA", "AC", "AT", "GC", "TG", "CA", "GA", "CT"]
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
        @test_throws ArgumentError GCATCodes._has_cycle_longer_than(g, 1)
    end


    @testset "invalid vertices but no edges" begin
        g = SimpleDiGraph(3)
        @test_throws ArgumentError GCATCodes._has_cycle_longer_than(g, 1)
    end


    @testset "invalid max_length argument" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        @test_throws ArgumentError GCATCodes._has_cycle_longer_than(g, -2)
    end
end


@testset "_is_codon_graphs_equal" begin
    @testset "happy path" begin
        codon_set_1 = LongDNA{4}.(["AAC", "ATG"])
        codon_set_2 = LongDNA{4}.(["AAC", "ATG"])
        cgd_1 = CodonGraphData(codon_set_1)
        cgd_2 = CodonGraphData(codon_set_2)
        @test GCATCodes._is_codon_graphs_equal(cgd_1, cgd_2) == true
    end


    @testset "sad path" begin
        codon_set_1 = LongDNA{4}.(["AGT", "GCA"])
        codon_set_2 = LongDNA{4}.(["AGT", "CGA"])
        cgd_1 = CodonGraphData(codon_set_1)
        cgd_2 = CodonGraphData(codon_set_2)
        @test GCATCodes._is_codon_graphs_equal(cgd_1, cgd_2) == false
    end


    @testset "same vertices but different edges" begin
        codon_set = LongDNA{4}.(["GAA", "GCA"])
        cgd_1 = CodonGraphData(codon_set)
        cgd_2 = CodonGraphData(codon_set)
        edge_label = ("A", "CA")
        GCATCodes._add_edge_by_label!(cgd_1, edge_label)
        @test GCATCodes._is_codon_graphs_equal(cgd_1, cgd_2) == false
    end


    @testset "same edges but different vertices" begin
        codon_set_1 = LongDNA{4}.(["AGT", "GCA"])
        codon_set_2 = LongDNA{4}.(["AGT", "GCA"])
        cgd_1 = CodonGraphData(codon_set_1)
        cgd_2 = CodonGraphData(codon_set_2)
        GCATCodes._add_vert_by_label!(cgd_1, "CT")
        @test GCATCodes._is_codon_graphs_equal(cgd_1, cgd_2) == false
    end
end


@testset "_has_edge_label" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CGT", "TTA"])
    cgd = CodonGraphData(codon_set)

    @testset "happy path" begin
        edge_label = ("C", "GT")
        @test GCATCodes._has_edge_label(cgd, edge_label) == true
    end


    @testset "sad path" begin
        edge_label = ("C", "TA")
        @test GCATCodes._has_edge_label(cgd, edge_label) == false
    end


    @testset "src_label not in vert_labels" begin
        edge_label = ("GA", "T")
        @test_throws ArgumentError GCATCodes._has_edge_label(cgd, edge_label)
    end


    @testset "dst_label not in vert_labels" begin
        edge_label = ("T", "AT")
        @test_throws ArgumentError GCATCodes._has_edge_label(cgd, edge_label)
    end


    @testset "src_label not in vert_idxs" begin
        edge_label = ("G", "AA")
        @test_throws ArgumentError GCATCodes._has_edge_label(cgd, edge_label)
    end


    @testset "dst_label not in vert_idxs" begin
        edge_label = ("T", "GG")
        @test_throws ArgumentError GCATCodes._has_edge_label(cgd, edge_label)
    end
end


@testset "_has_vert_label" begin
    # CodonGraphData object for all cases
    codon_set = LongDNA{4}.(["CTT", "ACT"])
    cgd = CodonGraphData(codon_set)


    @testset "happy path" begin
        vert_label = "CT"
        @test GCATCodes._has_vert_label(cgd, vert_label) == true
    end


    @testset "sad path" begin
        vert_label = "AG"
        @test GCATCodes._has_vert_label(cgd, vert_label) == false
    end
end