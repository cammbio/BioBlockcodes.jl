using GCATCodes
using Test


@testset "codon_set_to_str" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["CGT", "TCT"])
        @test codon_set_to_str(codon_set) == "\"CGT\", \"TCT\""
    end


    @testset "empty codon set" begin
        codon_set = LongDNA{4}[]
        @test_throws ArgumentError codon_set_to_str(codon_set)
    end


    @testset "invalid codon in codon set" begin
        codon_set = LongDNA{4}.(["AGT", "TGT", "GT"])
        @test_throws ArgumentError codon_set_to_str(codon_set)
    end
end


@testset "get_codon_set_from_line" begin
    @testset "happy path" begin
        line = "AAT|ACT,3|7"
        codon_set = get_codon_set_from_line(line)
        @test codon_set == LongDNA{4}.(["AAT", "ACT"])
    end


    @testset "empty line" begin
        line = ""
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "missing comma separator" begin
        line = "AAT|ACT"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "multiple comma separators" begin
        line = "AAT|ACT,3|7,extra"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "missing index list" begin
        line = "AAT|ACT,"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "non-numeric index list" begin
        line = "AAT|ACT,3|x"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "codon and index count mismatch" begin
        line = "AAT|ACT,3"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end


    @testset "index list does not match codon list" begin
        line = "AAT|ACT,4|7"
        @test_throws ArgumentError get_codon_set_from_line(line)
    end
end


@testset "_get_comb_from_codon_set" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["AAT", "ACT"])
        comb = GCATCodes._get_comb_from_codon_set(codon_set)
        @test comb == [3, 7]
    end


    @testset "duplicate codons" begin
        codon_set = LongDNA{4}.(["AAT", "AAT"])
        @test_throws ArgumentError GCATCodes._get_comb_from_codon_set(codon_set)
    end
end


@testset "_write_res" begin
    @testset "happy path" begin
        io = IOBuffer()
        codon_set = LongDNA{4}.(["AAT", "ACT"])
        comb = [3, 7]

        @test GCATCodes._write_res(io, codon_set, comb) == true
        @test String(take!(io)) == "AAT|ACT,3|7\n"
    end


    @testset "empty codon_set" begin
        io = IOBuffer()
        codon_set = LongDNA{4}[]
        comb = [1]
        @test_throws ArgumentError GCATCodes._write_res(io, codon_set, comb)
    end


    @testset "invalid comb" begin
        io = IOBuffer()
        codon_set = LongDNA{4}.(["AAT", "ACT"])
        comb = [7, 3]
        @test_throws ArgumentError GCATCodes._write_res(io, codon_set, comb)
    end
end
