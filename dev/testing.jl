using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs
using BioSequences
using NetworkLayout
using Base.Threads
using BenchmarkTools

# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
global_logger(ConsoleLogger(Logging.Info)) # deactivate

const ALL_CODONS =
    LongDNA{
        4,
    }.([
        "AAC",
        "AAG",
        "AAT",
        "ACA",
        "ACC",
        "ACG",
        "ACT",
        "AGA",
        "AGC",
        "AGG",
        "AGT",
        "ATA",
        "ATC",
        "ATG",
        "ATT",
        "CAA",
        "CAC",
        "CAG",
        "CAT",
        "CCA",
        "CCG",
        "CCT",
        "CGA",
        "CGC",
        "CGG",
        "CGT",
        "CTA",
        "CTC",
        "CTG",
        "CTT",
        "GAA",
        "GAC",
        "GAG",
        "GAT",
        "GCA",
        "GCC",
        "GCG",
        "GCT",
        "GGA",
        "GGC",
        "GGT",
        "GTA",
        "GTC",
        "GTG",
        "GTT",
        "TAA",
        "TAC",
        "TAG",
        "TAT",
        "TCA",
        "TCC",
        "TCG",
        "TCT",
        "TGA",
        "TGC",
        "TGG",
        "TGT",
        "TTA",
        "TTC",
        "TTG",
    ])

const stop_flag = Base.Threads.Atomic{Bool}(false)
codons = GCATCodes.ALL_CODONS

function sort_by_indices(infile::AbstractString, outfile::AbstractString)
    lines = readlines(infile)

    sorted = sort(lines; lt = (a, b) -> begin
        idxs_a = parse.(Int, split(last(split(a, ',')), '|'))
        idxs_b = parse.(Int, split(last(split(b, ',')), '|'))
        idxs_a < idxs_b
    end)

    open(outfile, "w") do io
        write(io, join(sorted, "\n"))
    end
end

for i in 2:2
    orig_file = "files/results/result_$i.csv"
    copy_file = "files/results/sorted_result_$i.csv"
    # sort_by_indices(orig_file, copy_file)
    for line in eachline(copy_file)
        combination = extract_csv_column(line, 2)
        println("Combination indices: ", combination)
    end
end

# check if there is an empty line in a file
function check_empty_lines_in_file(file_path::String)
    open(file_path, "r") do f
        line_number = 0
        for line in eachline(f)
            line_number += 1
            if isempty(strip(line))
                println("Empty line found at line number: $line_number")
            end
        end
    end
end

# function to get combination from codon set
function _get_combination_from_codon_set(codon_set::Vector{LongDNA{4}}, codons::Vector{LongDNA{4}})
    combination = Vector{Int}()
    codon_set_set = Set(codon_set)
    @inbounds for (i, c) in enumerate(codons)
        if c in codon_set_set
            push!(combination, i)
        end
    end
    return combination
end

open("files/results/test.txt", "r") do res
    counter = 0
    for line in eachline(res)
        if counter == 100_000
            break
        end
        codon_set = result_to_codon_set(line)
        # println("Codon set: $codon_set")
        counter += 1
    end
end

open("files/results/test.txt", "w") do in
    open("files/results/result_8.txt", "r") do out
        for line in eachline(out)
            println(in, line)
        end
    end
end

open("files/results/result_5.txt", "r") do res
    counter = 0
    false_positives = 0
    for line in eachline(res)
        # if counter == 10
        #     break
        # end
        counter += 1
        if counter % 100000 == 0
            println("counter: $counter")
        end

        codon_set = result_to_codon_set(line)
        data = CodonGraphData(codon_set)
        construct_graph_data!(data)
        _expand_graph(data)

        if get_max_cycle_length(data.graph; show_debug = false) > 2
            println("Max cycle length > 2 found!")
            println("codon set: $codon_set in line $counter is FALSE POSITIVE")
            false_positives += 1
        end
    end
    println("FINISHED, counter: $counter, false_positives: $false_positives")
end

function subsets(codon_set::Vector{LongDNA{4}})
    n = length(codon_set)
    result = Vector{Vector{LongDNA{4}}}()
    for i in 1:(2^n - 1)
        subset = Vector{LongDNA{4}}()
        for j in 1:n
            if (i >> (j - 1)) & 1 == 1
                push!(subset, codon_set[j])
            end
        end
        push!(result, subset)
    end
    return result
end

# function to get processed count from current combination
function _get_processed_count_from_combination(comb::Vector{Int}; n::Int = 60)
    k = length(comb)
    rank = 0
    prev = 0
    @inbounds for (i, ci) in enumerate(comb)
        for x in (prev + 1):(ci - 1)
            rank += binomial(n - x, k - i)
        end
        prev = ci
    end
    return rank
end
