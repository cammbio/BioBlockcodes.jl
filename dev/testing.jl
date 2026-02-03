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

const stop_flag = Base.Threads.Atomic{Bool}(false)
codons = GCATCodes.ALL_CODONS











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

check_empty_lines_in_file("files/results/result_7.txt")
check_empty_lines_in_file("files/results/result_8.txt")
check_empty_lines_in_file("files/results/result_9.txt")
check_empty_lines_in_file("files/results/result_10.txt")
check_empty_lines_in_file("files/results/result_11.txt")
check_empty_lines_in_file("files/results/result_12.txt")
check_empty_lines_in_file("files/results/result_13.txt")
check_empty_lines_in_file("files/results/result_14.txt")
check_empty_lines_in_file("files/results/result_15.txt")
check_empty_lines_in_file("files/results/result_16.txt")
check_empty_lines_in_file("files/results/result_17.txt")
check_empty_lines_in_file("files/results/result_18.txt")
check_empty_lines_in_file("files/results/result_19.txt")
check_empty_lines_in_file("files/results/result_20.txt")

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
# -------------------------------------------------- BENCHMARKING --------------------------------------------------
const MICRO_ITERS = 1_000_000
codons = GCATCodes.ALL_CODONS
rotation_masks = GCATCodes._build_rotation_masks(codons)

index_map = Dict{LongDNA{4}, Int}()
@inbounds for (i, c) in enumerate(codons)
    index_map[c] = i
end

base_idx = 1
base = codons[base_idx]
idx1 = index_map[left_shift_codon(base, 1)]
idx2 = index_map[left_shift_codon(base, 2)]
extra_idx = base_idx == 1 ? 2 : 1
if extra_idx == idx1 || extra_idx == idx2
    extra_idx = 3
end
combination = [base_idx, idx1, idx2, extra_idx]

function _reject_plain_once(codons, combination)
    codon_set = codons[combination]
    return _contains_codon_rotation(codon_set)
end

function _reject_mask_once(rotation_masks, combination)
    comb_mask = GCATCodes._combination_to_mask(combination)
    return GCATCodes._mask_contains_codon_rotation(comb_mask, combination, rotation_masks)
end

function _micro_reject_plain(codons, combination, iters)
    hits = 0
    for _ in 1:iters
        _reject_plain_once(codons, combination) && (hits += 1)
    end
    return hits
end

function _micro_reject_mask(rotation_masks, combination, iters)
    hits = 0
    for _ in 1:iters
        _reject_mask_once(rotation_masks, combination) && (hits += 1)
    end
    return hits
end

trial_plain = @benchmark _micro_reject_plain($codons, $combination, $MICRO_ITERS) samples = 500 evals = 1
trial_mask =
    @benchmark _micro_reject_mask($rotation_masks, $combination, $MICRO_ITERS) samples = 500 evals = 1
m_plain = median(trial_plain).time
m_mask = median(trial_mask).time
println("plain: $(m_plain) ns for $MICRO_ITERS iterations")
println("mask : $(m_mask) ns for $MICRO_ITERS iterations")
if m_mask < m_plain
    println("speedup: $(round(m_plain / m_mask, digits = 2))x (mask faster)")
else
    println("speedup: $(round(m_mask / m_plain, digits = 2))x (plain faster)")
end
