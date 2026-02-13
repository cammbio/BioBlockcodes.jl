using Test
using GCATCodes


println("Running tests for InputOutputUtilities.jl...")


@testset "codon_set_to_str" begin
    # happy path
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["CGT", "TCT"])
        @test codon_set_to_str(codon_set) == "CGT,TCT"
    end


    # empty codon set
    @testset "empty codon set" begin
        codon_set = LongDNA{4}[]
        @test codon_set_to_str(codon_set) == ""
    end
end


# get_codon_set_from_line
@testset "get_codon_set_from_line" begin
    # happy path
    @testset "happy path" begin
        line = "AAT,ACT|4,5"
        codon_set = get_codon_set_from_line(line)
        @test codon_set == LongDNA{4}.(["AAT", "ACT"])
    end
end