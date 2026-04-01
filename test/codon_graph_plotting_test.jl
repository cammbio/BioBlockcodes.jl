using BioSequences
using CairoMakie
using BioBlockcodes
using Graphs
using Test


@testset "plot_codon_graph" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        fig = plot_codon_graph(cgd)
        @test isa(fig, Figure)
    end


    @testset "custom fig_size" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        fig_size = (1200, 700)
        fig = plot_codon_graph(cgd; fig_size = fig_size)
        @test Tuple(fig.scene.viewport[].widths) == fig_size
    end


    @testset "sets graph title on axis" begin
        graph_title = "My Test Graph"
        codon_set = LongDNA{4}.(["GAC", "TTA"])
        cgd = CodonGraphData(codon_set; graph_title = graph_title)
        fig = plot_codon_graph(cgd)
        # find the axis in the figure content and check its title
        axis = only(filter(obj -> obj isa Axis, fig.content))
        @test axis.title[] == graph_title
    end


    @testset "empty title accepted" begin
        codon_set = LongDNA{4}.(["TCC", "GTA"])
        cgd = CodonGraphData(codon_set; graph_title = "")
        fig = plot_codon_graph(cgd)
        @test isa(fig, Figure)
    end


    @testset "does not mutate input cgd" begin
        codon_set = LongDNA{4}.(["AAC", "TGT"])
        cgd = CodonGraphData(codon_set; graph_title = "Immutable Input Test")

        orig_codon_set = copy(cgd.codon_set)
        orig_edge_labels = copy(cgd.edge_labels)
        orig_vert_labels = copy(cgd.vert_labels)
        orig_vert_idxs = copy(cgd.vert_idxs)
        orig_graph_title = cgd.graph_title
        orig_nv = nv(cgd.graph)
        orig_ne = ne(cgd.graph)
        orig_edges = collect(edges(cgd.graph))

        plot_codon_graph(cgd)

        @test cgd.codon_set == orig_codon_set
        @test cgd.edge_labels == orig_edge_labels
        @test cgd.vert_labels == orig_vert_labels
        @test cgd.vert_idxs == orig_vert_idxs
        @test cgd.graph_title == orig_graph_title
        @test nv(cgd.graph) == orig_nv
        @test ne(cgd.graph) == orig_ne
        @test collect(edges(cgd.graph)) == orig_edges
    end
end


@testset "plot_multiple_codon_graphs" begin
    @testset "happy path" begin
        codon_set_1 = LongDNA{4}.(["ATG", "CAT"])
        codon_set_2 = LongDNA{4}.(["GAC", "TTA"])
        cgd_1 = CodonGraphData(codon_set_1; graph_title = "Graph 1")
        cgd_2 = CodonGraphData(codon_set_2; graph_title = "Graph 2")
        cgd_list = [cgd_1, cgd_2]
        fig = plot_multiple_codon_graphs(cgd_list)
        @test isa(fig, Figure)
    end


    @testset "custom fig_size" begin
        codon_set_1 = LongDNA{4}.(["ATG", "CAT"])
        codon_set_2 = LongDNA{4}.(["GAC", "TTA"])
        cgd_1 = CodonGraphData(codon_set_1; graph_title = "Graph 1")
        cgd_2 = CodonGraphData(codon_set_2; graph_title = "Graph 2")
        cgd_list = [cgd_1, cgd_2]
        fig_size = (1400, 800)
        fig = plot_multiple_codon_graphs(cgd_list; fig_size = fig_size)
        @test Tuple(fig.scene.viewport[].widths) == fig_size
    end


    @testset "creates one axis per graph with matching titles" begin
        codon_set_1 = LongDNA{4}.(["ATG", "CAT"])
        codon_set_2 = LongDNA{4}.(["TCC", "GTA"])
        codon_set_3 = LongDNA{4}.(["AAC", "TGT"])
        cgd_1 = CodonGraphData(codon_set_1; graph_title = "Graph A")
        cgd_2 = CodonGraphData(codon_set_2; graph_title = "Graph B")
        cgd_3 = CodonGraphData(codon_set_3; graph_title = "Graph C")
        cgd_list = [cgd_1, cgd_2, cgd_3]
        fig = plot_multiple_codon_graphs(cgd_list)
        # find all axes in the figure content and check their titles
        axes = filter(obj -> obj isa Axis, fig.content)
        @test length(axes) == length(cgd_list)
        @test [ax.title[] for ax in axes] == getfield.(cgd_list, :graph_title)
    end


    @testset "fig title label is optional" begin
        codon_set_1 = LongDNA{4}.(["ATG", "CAT"])
        codon_set_2 = LongDNA{4}.(["GAC", "TTA"])
        cgd_1 = CodonGraphData(codon_set_1; graph_title = "Graph 1")
        cgd_2 = CodonGraphData(codon_set_2; graph_title = "Graph 2")
        cgd_list = [cgd_1, cgd_2]

        fig_without_title = plot_multiple_codon_graphs(cgd_list; fig_title = nothing)
        labels_without_title = filter(obj -> obj isa Label, fig_without_title.content)
        @test length(labels_without_title) == 0

        fig_with_empty_title = plot_multiple_codon_graphs(cgd_list; fig_title = "")
        labels_with_empty_title = filter(obj -> obj isa Label, fig_with_empty_title.content)
        @test length(labels_with_empty_title) == 0

        fig_with_title = plot_multiple_codon_graphs(cgd_list; fig_title = "Overview")
        labels_with_title = filter(obj -> obj isa Label, fig_with_title.content)
        @test length(labels_with_title) == 1
        @test only(labels_with_title).text[] == "Overview"
    end


    @testset "empty cgd_list throws ArgumentError" begin
        @test_throws ArgumentError plot_multiple_codon_graphs(CodonGraphData[])
    end


    @testset "does not mutate input cgd_list and contained cgd objects" begin
        codon_set_1 = LongDNA{4}.(["ATG", "CAT"])
        codon_set_2 = LongDNA{4}.(["GAC", "TTA"])
        cgd_1 = CodonGraphData(codon_set_1; graph_title = "Graph 1")
        cgd_2 = CodonGraphData(codon_set_2; graph_title = "Graph 2")
        cgd_list = [cgd_1, cgd_2]

        orig_cgd_list = copy(cgd_list)
        orig_state = [
            (
                codon_set = copy(d.codon_set),
                edge_labels = copy(d.edge_labels),
                vert_labels = copy(d.vert_labels),
                vert_idxs = copy(d.vert_idxs),
                graph_title = d.graph_title,
                nv = nv(d.graph),
                ne = ne(d.graph),
                edges = collect(edges(d.graph)),
            ) for d in cgd_list
        ]

        plot_multiple_codon_graphs(cgd_list; fig_title = "Overview")

        @test cgd_list == orig_cgd_list
        for (idx, cgd) in enumerate(cgd_list)
            @test cgd.codon_set == orig_state[idx].codon_set
            @test cgd.edge_labels == orig_state[idx].edge_labels
            @test cgd.vert_labels == orig_state[idx].vert_labels
            @test cgd.vert_idxs == orig_state[idx].vert_idxs
            @test cgd.graph_title == orig_state[idx].graph_title
            @test nv(cgd.graph) == orig_state[idx].nv
            @test ne(cgd.graph) == orig_state[idx].ne
            @test collect(edges(cgd.graph)) == orig_state[idx].edges
        end
    end
end


@testset "_get_col_count" begin
    @testset "normal values" begin
        @test BioBlockcodes._get_col_count(1) == 1
        @test BioBlockcodes._get_col_count(2) == 2
        @test BioBlockcodes._get_col_count(3) == 2
        @test BioBlockcodes._get_col_count(5) == 3
    end


    @testset "perfect squares" begin
        @test BioBlockcodes._get_col_count(4) == 2
        @test BioBlockcodes._get_col_count(9) == 3
        @test BioBlockcodes._get_col_count(16) == 4
    end


    @testset "zero graph count" begin
        @test_throws ArgumentError BioBlockcodes._get_col_count(0)
    end


    @testset "negative graph count" begin
        @test_throws ArgumentError BioBlockcodes._get_col_count(-1)
    end
end
