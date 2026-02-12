using BioSequences
using Graphs
using Test

println("Running GraphUtils tests...")


@testset "construct_graph_data!" begin
    # case 1 with accepted codon set with optional title
    @testset "case 1" begin
        codon_set = LongDNA{4}.(["CCA", "ATG"])
        plot_title = "Test title"
        data = CodonGraphData(codon_set; plot_title = plot_title)
        construct_graph_data!(data; show_debug = false)

        # test graph content
        # check that right amount of vertices and edges are in the graph
        @test nv(data.graph) == 7
        @test ne(data.graph) == 4
        @test length(data.all_vertex_labels) == nv(data.graph)
        @test length(data.vertex_index) == nv(data.graph)
        @test length(data.all_edge_labels) == ne(data.graph)
        # check that all vertices are in the graph
        for vertex in data.all_vertex_labels
            @test haskey(data.vertex_index, vertex)
            @test has_vertex(data.graph, data.vertex_index[vertex])
        end
        # check that all edges are in the graph
        for (src_label, dst_label) in data.all_edge_labels
            @test haskey(data.vertex_index, src_label)
            @test haskey(data.vertex_index, dst_label)
            @test has_edge(data.graph, data.vertex_index[src_label], data.vertex_index[dst_label])
        end

        # test fields of CodonGraphData
        @test data.codon_set == codon_set
        @test data.all_vertex_labels == ["A", "C", "G", "AT", "CA", "CC", "TG"]
        @test data.base_vertex_labels == data.all_vertex_labels
        @test isempty(data.added_vertex_labels)
        @test data.all_edge_labels == [("C", "CA"), ("CC", "A"), ("A", "TG"), ("AT", "G")]
        @test data.base_edge_labels == data.all_edge_labels
        @test isempty(data.added_edge_labels)
        @test data.vertex_index ==
              Dict("AT" => 4, "A" => 1, "CA" => 5, "C" => 2, "TG" => 7, "CC" => 6, "G" => 3)
        @test data.base_vertex_labels[data.vertex_index["CA"]] == "CA"
        @test data.plot_title == plot_title
    end

    # case 2 with accepted codon set without optional title
    @testset "case 2" begin
        codon_set = LongDNA{4}.(["GTA", "GTT", "GCA"])
        data = CodonGraphData(codon_set)
        construct_graph_data!(data; show_debug = false)

        # test graph content
        # check that right amount of vertices and edges are in the graph
        @test nv(data.graph) == 8
        @test ne(data.graph) == 6
        @test length(data.all_vertex_labels) == nv(data.graph)
        @test length(data.vertex_index) == nv(data.graph)
        @test length(data.all_edge_labels) == ne(data.graph)
        # check that all vertices are in the graph
        for vertex in data.all_vertex_labels
            @test haskey(data.vertex_index, vertex)
            @test has_vertex(data.graph, data.vertex_index[vertex])
        end
        # check that all edges are in the graph
        for (src_label, dst_label) in data.all_edge_labels
            @test haskey(data.vertex_index, src_label)
            @test haskey(data.vertex_index, dst_label)
            @test has_edge(data.graph, data.vertex_index[src_label], data.vertex_index[dst_label])
        end

        # test fields of CodonGraphData
        @test data.codon_set == codon_set
        @test data.all_vertex_labels == ["A", "G", "T", "CA", "GC", "GT", "TA", "TT"]
        @test data.base_vertex_labels == data.all_vertex_labels
        @test isempty(data.added_vertex_labels)
        @test data.all_edge_labels ==
              [("G", "TA"), ("GT", "A"), ("G", "TT"), ("GT", "T"), ("G", "CA"), ("GC", "A")]
        @test data.base_edge_labels == data.all_edge_labels
        @test isempty(data.added_edge_labels)
        @test data.vertex_index ==
              Dict("A" => 1, "T" => 3, "CA" => 4, "GC" => 5, "G" => 2, "GT" => 6, "TA" => 7, "TT" => 8)
        @test data.base_vertex_labels[data.vertex_index["GC"]] == "GC"
        @test data.plot_title == ""
    end

end


@testset "add_vertices_by_codon_set!" begin
    # case 1
    @testset "case 1" begin
        codon_set = LongDNA{4}.(["AAC", "GTT"])
        data = CodonGraphData(codon_set)
        add_vertices_by_codon_set!(data; show_debug = false)

        # test graph content
        @test nv(data.graph) == length(data.all_vertex_labels)
        @test ne(data.graph) == 0

        # test fields of CodonGraphData
        @test data.codon_set == codon_set
        @test data.all_vertex_labels == ["A", "C", "G", "T", "AA", "AC", "GT", "TT"]
        @test data.base_vertex_labels == data.all_vertex_labels
        @test isempty(data.added_vertex_labels)
        @test isempty(data.all_edge_labels)
        @test isempty(data.base_edge_labels)
        @test isempty(data.added_edge_labels)
        @test isempty(data.vertex_index)

        # test sorting
        @test issorted(data.all_vertex_labels, by = x -> (length(x), x))
    end

    # case 2 with overlapping codons
    @testset "case 2" begin
        codon_set = LongDNA{4}.(["AAC", "CAA"])
        data = CodonGraphData(codon_set)
        add_vertices_by_codon_set!(data; show_debug = false)

        # test graph content
        @test nv(data.graph) == length(data.all_vertex_labels)
        @test ne(data.graph) == 0

        # test fields of CodonGraphData
        @test data.codon_set == codon_set
        @test data.all_vertex_labels == ["A", "C", "AA", "AC", "CA"]
        @test data.base_vertex_labels == data.all_vertex_labels
        @test isempty(data.added_vertex_labels)
        @test isempty(data.all_edge_labels)
        @test isempty(data.base_edge_labels)
        @test isempty(data.added_edge_labels)
        @test isempty(data.vertex_index)

        # test sorting
        @test issorted(data.all_vertex_labels, by = x -> (length(x), x))
    end
end



@testset "add_edges_by_codon_set!" begin
    codon_set = LongDNA{4}.(["ATG", "TTC"])
    data = CodonGraphData(codon_set)
    add_vertices_by_codon_set!(data; show_debug = false)
    data.vertex_index = Dict(label => index for (index, label) in enumerate(data.all_vertex_labels))
    add_edges_by_codon_set!(data; show_debug = false)

    # test graph content
    @test nv(data.graph) == length(data.all_vertex_labels)
    @test ne(data.graph) == 2 * length(codon_set)
    @test ne(data.graph) == length(data.all_edge_labels)

    # test fields of CodonGraphData
    @test data.codon_set == codon_set
    @test data.all_vertex_labels == ["A", "C", "G", "T", "AT", "TC", "TG", "TT"]
    @test data.base_vertex_labels == data.all_vertex_labels
    @test isempty(data.added_vertex_labels)
    @test data.all_edge_labels == [("A", "TG"), ("AT", "G"), ("T", "TC"), ("TT", "C")]
    @test data.base_edge_labels == data.all_edge_labels
    @test isempty(data.added_edge_labels)
    @test data.vertex_index ==
          Dict("AT" => 5, "A" => 1, "T" => 4, "C" => 2, "TG" => 7, "TC" => 6, "G" => 3, "TT" => 8)
    @test data.base_vertex_labels[data.vertex_index["TG"]] == "TG"
end


