using BioSequences
using GCATCodes
using Test


@testset "get_comp_rev_codon_set" begin
    @testset "happy path single codon" begin
        codon_set = LongDNA{4}.(["GCT"])
        comp_rev_codon_set = get_comp_rev_codon_set(codon_set)
        @test comp_rev_codon_set == LongDNA{4}.(["AGC"])
    end


    @testset "happy path multiple codons" begin
        codon_set = LongDNA{4}.(["CGA", "TTA", "GGC"])
        comp_rev_codon_set = get_comp_rev_codon_set(codon_set)
        @test comp_rev_codon_set == LongDNA{4}.(["TCG", "TAA", "GCC"])
    end


    @testset "returns same length as input" begin
        codon_set = LongDNA{4}.(["ACC", "TGA", "CGT", "AAT"])
        comp_rev_codon_set = get_comp_rev_codon_set(codon_set)
        @test length(comp_rev_codon_set) == length(codon_set)
    end


    @testset "does not mutate input codon_set" begin
        codon_set = LongDNA{4}.(["GTT", "CCA", "AGC"])
        copy_codon_set = copy(codon_set)
        comp_rev_codon_set = get_comp_rev_codon_set(codon_set)
        @test codon_set == copy_codon_set
    end


    @testset "involution property" begin
        codon_set = LongDNA{4}.(["TGC", "CAA", "GAT", "CTT"])
        comp_rev_codon_set = get_comp_rev_codon_set(codon_set)
        comp_rev_comp_rev_codon_set = get_comp_rev_codon_set(comp_rev_codon_set)
        @test comp_rev_comp_rev_codon_set == codon_set
    end
end


@testset "left_shift_codon" begin
    @testset "happy path shift by 1" begin
        codon = LongDNA{4}("ATG")
        shifted_codon = left_shift_codon(codon, 1)
        @test shifted_codon == LongDNA{4}("TGA")
    end


    @testset "happy path shift by 2" begin
        codon = LongDNA{4}("CGA")
        shifted_codon = left_shift_codon(codon, 2)
        @test shifted_codon == LongDNA{4}("ACG")
    end


    @testset "happy path shift by 0" begin
        codon = LongDNA{4}("GCT")
        shifted_codon = left_shift_codon(codon, 0)
        @test shifted_codon == codon
        @test shifted_codon !== codon
    end


    @testset "happy path shift by codon length" begin
        codon = LongDNA{4}("TGA")
        shifted_codon = left_shift_codon(codon, 3)
        @test shifted_codon == codon
        @test shifted_codon !== codon
    end


    @testset "shift wraps with modulo codon length" begin
        codon = LongDNA{4}("TTC")
        shift_by_1 = left_shift_codon(codon, 1)
        shift_by_4 = left_shift_codon(codon, 4)
        @test shift_by_4 == shift_by_1
    end


    @testset "negative shift throws ArgumentError" begin
        codon = LongDNA{4}("AAT")
        @test_throws ArgumentError left_shift_codon(codon, -1)
    end


    @testset "does not mutate input codon" begin
        codon = LongDNA{4}("CGC")
        copy_codon = copy(codon)
        shifted_codon = left_shift_codon(codon, 2)
        @test codon == copy_codon
    end
end


@testset "left_shift_codon_set" begin
    @testset "happy path shift by 1" begin
        codon_set = LongDNA{4}.(["ATG", "CCA", "GTT"])
        shifted_codon_set = left_shift_codon_set(codon_set, 1)
        @test shifted_codon_set == LongDNA{4}.(["TGA", "CAC", "TTG"])
    end


    @testset "happy path shift by 2" begin
        codon_set = LongDNA{4}.(["AGC", "TTA", "CGA"])
        shifted_codon_set = left_shift_codon_set(codon_set, 2)
        @test shifted_codon_set == LongDNA{4}.(["CAG", "ATT", "ACG"])
    end


    @testset "happy path shift by 0" begin
        codon_set = LongDNA{4}.(["GCT", "AAT", "TGC"])
        shifted_codon_set = left_shift_codon_set(codon_set, 0)
        @test shifted_codon_set == codon_set
        @test shifted_codon_set !== codon_set
    end


    @testset "shift wraps with modulo codon length" begin
        codon_set = LongDNA{4}.(["ATG", "CAA", "TTC"])
        shift_by_1 = left_shift_codon_set(codon_set, 1)
        shift_by_4 = left_shift_codon_set(codon_set, 4)
        @test shift_by_4 == shift_by_1
    end


    @testset "negative shift throws ArgumentError" begin
        codon_set = LongDNA{4}.(["ATG", "CAA", "TTC"])
        @test_throws ArgumentError left_shift_codon_set(codon_set, -1)
    end


    @testset "does not mutate input codon_set" begin
        codon_set = LongDNA{4}.(["TGA", "CGC", "AAT"])
        copy_codon_set = copy(codon_set)
        shifted_codon_set = left_shift_codon_set(codon_set, 2)
        @test codon_set == copy_codon_set
    end
end


@testset "_get_comp_base" begin
    @testset "happy path" begin
        @test GCATCodes._get_comp_base(DNA_A) == DNA_T
        @test GCATCodes._get_comp_base(DNA_C) == DNA_G
        @test GCATCodes._get_comp_base(DNA_G) == DNA_C
        @test GCATCodes._get_comp_base(DNA_T) == DNA_A
    end


    @testset "invalid base throws ArgumentError" begin
        @test_throws ArgumentError GCATCodes._get_comp_base(DNA_N)
    end
end


@testset "_get_comp_codon" begin
    @testset "happy path" begin
        codon = LongDNA{4}("ATG")
        comp_codon = GCATCodes._get_comp_codon(codon)
        @test comp_codon == LongDNA{4}("TAC")
    end


    @testset "does not mutate input codon" begin
        codon = LongDNA{4}("CGA")
        copy_codon = copy(codon)
        comp_codon = GCATCodes._get_comp_codon(codon)
        @test codon == copy_codon
    end
end


@testset "_get_rev_codon" begin
    @testset "happy path" begin
        codon = LongDNA{4}("ATG")
        rev_codon = GCATCodes._get_rev_codon(codon)
        @test rev_codon == LongDNA{4}("GTA")
    end


    @testset "does not mutate input codon" begin
        codon = LongDNA{4}("TTC")
        copy_codon = copy(codon)
        rev_codon = GCATCodes._get_rev_codon(codon)
        @test codon == copy_codon
    end
end


@testset "_get_comp_codon_set" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CGA", "TTC"])
        comp_codon_set = GCATCodes._get_comp_codon_set(codon_set)
        @test comp_codon_set == LongDNA{4}.(["TAC", "GCT", "AAG"])
    end


    @testset "returns same length as input" begin
        codon_set = LongDNA{4}.(["GCT", "AAT", "TGC", "CCA"])
        comp_codon_set = GCATCodes._get_comp_codon_set(codon_set)
        @test length(comp_codon_set) == length(codon_set)
    end


    @testset "does not mutate input codon_set" begin
        codon_set = LongDNA{4}.(["AGC", "TTA", "CGA"])
        copy_codon_set = copy(codon_set)
        comp_codon_set = GCATCodes._get_comp_codon_set(codon_set)
        @test codon_set == copy_codon_set
    end
end


@testset "_get_rev_codon_set" begin
    @testset "happy path" begin
        codon_set = LongDNA{4}.(["ATG", "CGA", "TTC"])
        rev_codon_set = GCATCodes._get_rev_codon_set(codon_set)
        @test rev_codon_set == LongDNA{4}.(["GTA", "AGC", "CTT"])
    end


    @testset "returns same length as input" begin
        codon_set = LongDNA{4}.(["GCT", "AAT", "TGC", "CCA"])
        rev_codon_set = GCATCodes._get_rev_codon_set(codon_set)
        @test length(rev_codon_set) == length(codon_set)
    end


    @testset "does not mutate input codon_set" begin
        codon_set = LongDNA{4}.(["AGC", "TTA", "CGA"])
        copy_codon_set = copy(codon_set)
        rev_codon_set = GCATCodes._get_rev_codon_set(codon_set)
        @test codon_set == copy_codon_set
    end
end
