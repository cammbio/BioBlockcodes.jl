using BioSequences
using BioBlockcodes
using Graphs
using Test


@testset "_validate_cgd" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        @test_nowarn BioBlockcodes._validate_cgd(cgd)
    end


    @testset "invalid edges" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        empty!(cgd.edge_labels)
        @test_throws ArgumentError BioBlockcodes._validate_cgd(cgd)
    end


    @testset "invalid vertices" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        empty!(cgd.vert_labels)
        @test_throws ArgumentError BioBlockcodes._validate_cgd(cgd)
    end
end


@testset "_validate_ckp_file" begin
    @testset "happy path" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "ckp.csv")
            write(ckp_path, "comb_size,2\nnext_line,15\nstatus,unfinished\n")
            @test_nowarn BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "file not found" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "missing.csv")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "empty file" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "empty.csv")
            touch(ckp_path)
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid line count" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nnext_line,15\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid line format" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2,extra\nnext_line,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid whitespace after comma" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size, 2\nnext_line,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid whitespace before comma" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size ,2\nnext_line,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid key" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nline,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "duplicate key" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\ncomb_size,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "non-integer comb_size" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,two\nnext_line,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "empty comb_size value" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,\nnext_line,15\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "non-integer next_line" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nnext_line,fifteen\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "empty next_line value" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nnext_line,\nstatus,unfinished\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "invalid status value" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nnext_line,15\nstatus,done\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


    @testset "empty status value" begin
        mktempdir() do temp_dir
            ckp_path = joinpath(temp_dir, "bad.csv")
            write(ckp_path, "comb_size,2\nnext_line,15\nstatus,\n")
            @test_throws ArgumentError BioBlockcodes._validate_ckp_file(ckp_path)
        end
    end


end


@testset "_validate_ckp" begin
    @testset "happy path unfinished" begin
        ckp = (comb_size = 2, next_line = 15, status = "unfinished")
        @test_nowarn BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "happy path finished" begin
        ckp = (comb_size = 2, next_line = 0, status = "finished")
        @test_nowarn BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "comb size mismatch" begin
        ckp = (comb_size = 3, next_line = 15, status = "unfinished")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "current line below 1" begin
        ckp = (comb_size = 2, next_line = 0, status = "unfinished")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "current line above max lines" begin
        ckp = (comb_size = 2, next_line = 61, status = "unfinished")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "invalid status" begin
        ckp = (comb_size = 2, next_line = 15, status = "done")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "finished status with non-max line" begin
        ckp = (comb_size = 2, next_line = 59, status = "finished")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "finished status with max line" begin
        ckp = (comb_size = 2, next_line = 60, status = "finished")
        @test_throws ArgumentError BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end


    @testset "unfinished status with max line" begin
        ckp = (comb_size = 2, next_line = 60, status = "unfinished")
        @test_nowarn BioBlockcodes._validate_ckp(ckp, "files/tests/checkpoints/ckp_2.csv", 2, 60)
    end
end


@testset "_validate_comb" begin
    @testset "happy path" begin
        comb = [1, 5, 12]
        @test_nowarn BioBlockcodes._validate_comb(comb)
    end


    @testset "empty comb" begin
        comb = Int[]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end


    @testset "invalid index" begin
        comb = [1, 61]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end


    @testset "invalid zero index" begin
        comb = [0, 5, 12]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end


    @testset "invalid negative index" begin
        comb = [-2, 5, 12]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end


    @testset "duplicate indices" begin
        comb = [2, 2, 5]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end


    @testset "unsorted comb" begin
        comb = [1, 7, 3]
        @test_throws ArgumentError BioBlockcodes._validate_comb(comb)
    end
end


@testset "_validate_dir" begin
    @testset "happy path" begin
        mktempdir() do temp_dir
            file_path = joinpath(temp_dir, "out.csv")
            @test_nowarn BioBlockcodes._validate_dir(file_path)
        end
    end


    @testset "directory not found" begin
        temp_dir = mktempdir()
        file_path = joinpath(temp_dir, "missing", "out.csv")
        @test_throws ArgumentError BioBlockcodes._validate_dir(file_path)
    end
end


@testset "_validate_graph" begin
    @testset "happy path" begin
        graph = SimpleDiGraph(2)
        add_edge!(graph, 1, 2)
        @test_nowarn BioBlockcodes._validate_graph(graph)
    end


    @testset "no edges" begin
        graph = SimpleDiGraph(3)
        @test_throws ArgumentError BioBlockcodes._validate_graph(graph)
    end


    @testset "no vertices" begin
        graph = SimpleDiGraph(0)
        @test_throws ArgumentError BioBlockcodes._validate_graph(graph)
    end
end


@testset "_validate_label" begin
    @testset "happy path" begin
        @test_nowarn BioBlockcodes._validate_label("A")
        @test_nowarn BioBlockcodes._validate_label("TG")
    end


    @testset "empty label" begin
        @test_throws ArgumentError BioBlockcodes._validate_label("")
    end


    @testset "invalid length" begin
        @test_throws ArgumentError BioBlockcodes._validate_label("ATG")
    end


    @testset "invalid character" begin
        @test_throws ArgumentError BioBlockcodes._validate_label("AN")
    end
end


@testset "_validate_codon" begin
    @testset "happy path" begin
        codon = LongDNA{4}("ATG")
        @test_nowarn BioBlockcodes._validate_codon(codon)
    end


    @testset "empty codon" begin
        codon = LongDNA{4}("")
        @test_throws ArgumentError BioBlockcodes._validate_codon(codon)
    end


    @testset "invalid codon" begin
        codon = LongDNA{4}("AA")
        @test_throws ArgumentError BioBlockcodes._validate_codon(codon)
    end
end


@testset "_validate_codon_set" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        @test_nowarn BioBlockcodes._validate_codon_set(codon_set)
    end


    @testset "empty codon_set" begin
        codon_set = LongDNA{4}[]
        @test_throws ArgumentError BioBlockcodes._validate_codon_set(codon_set)
    end


    @testset "duplicate codons" begin
        codon_set = LongDNA{4}.(["ATG", "ATG"])
        @test_throws ArgumentError BioBlockcodes._validate_codon_set(codon_set)
    end


    @testset "invalid codon in set" begin
        codon_set = LongDNA{4}.(["ATG", "AA"])
        @test_throws ArgumentError BioBlockcodes._validate_codon_set(codon_set)
    end
end


@testset "_validate_edges" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        @test_nowarn BioBlockcodes._validate_edges(cgd)
    end


    @testset "empty edge_labels" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        empty!(cgd.edge_labels)
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "edge_labels and graph edge count mismatch" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        pop!(cgd.edge_labels)
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "graph has no edges" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        cgd.graph = SimpleDiGraph(nv(cgd.graph))
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "duplicate edge labels" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        push!(cgd.edge_labels, first(cgd.edge_labels))
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "edge label source missing in vert_idxs" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        dst_label = cgd.edge_labels[1][2]
        cgd.edge_labels[1] = ("AA", dst_label)
        delete!(cgd.vert_idxs, "AA")
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "edge label source missing in vert_labels" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        dst_label = cgd.edge_labels[1][2]
        cgd.edge_labels[1] = ("AA", dst_label)
        cgd.vert_idxs["AA"] = 1
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end


    @testset "edge label not present in graph" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        labels = cgd.vert_labels
        bad_edge_label = nothing
        for src_label in labels
            for dst_label in labels
                if !BioBlockcodes._has_edge_label(cgd, (src_label, dst_label))
                    bad_edge_label = (src_label, dst_label)
                    break
                end
            end
            bad_edge_label !== nothing && break
        end
        cgd.edge_labels[1] = bad_edge_label
        @test_throws ArgumentError BioBlockcodes._validate_edges(cgd)
    end
end


@testset "_validate_vertices" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        @test_nowarn BioBlockcodes._validate_vertices(cgd)
    end


    @testset "empty vert_labels" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        empty!(cgd.vert_labels)
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "empty vert_idxs" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        empty!(cgd.vert_idxs)
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "graph has no vertices" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        cgd.graph = SimpleDiGraph(0)
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "duplicate vert_labels" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        push!(cgd.vert_labels, first(cgd.vert_labels))
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "invalid vert_label length" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        cgd.vert_labels[1] = "ATG"
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "invalid vert_label character" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        cgd.vert_labels[1] = "N"
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "vert_labels and graph vertex count mismatch" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        pop!(cgd.vert_labels)
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "vert_idxs and graph vertex count mismatch" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        delete!(cgd.vert_idxs, first(keys(cgd.vert_idxs)))
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "vert_idxs contains out of range index" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        label = first(keys(cgd.vert_idxs))
        cgd.vert_idxs[label] = nv(cgd.graph) + 1
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end


    @testset "vert_idxs and vert_labels mapping mismatch" begin
        codon_set = LongDNA{4}.(["ATG", "CAT"])
        cgd = CodonGraphData(codon_set)
        label = cgd.vert_labels[1]
        cgd.vert_idxs[label] = 2
        @test_throws ArgumentError BioBlockcodes._validate_vertices(cgd)
    end
end
