using Test
using GCATCodes


println("Running tests for InputOutputUtilities.jl...")


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
end