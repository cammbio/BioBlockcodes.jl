using Test
using BioSequences
using GCATCodes


println("Running tests for CodonGraphAnalysis.jl...")


# helper to create a codon graph for testing
function _create_data(codon_set::Vector{LongDNA{4}})
    return CodonGraphData(codon_set)
end


@testset ""
