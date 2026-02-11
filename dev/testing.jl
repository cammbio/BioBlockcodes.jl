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
# global_logger(ConsoleLogger(Logging.Info)) # deactivate

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

function grow_codon_results_all(src_codon_set::Vector{LongDNA{4}}; line_count::Int = 0, debug::Bool = false)
    data_list = Vector{CodonGraphData}()
    for i in 1:length(src_codon_set)
        codon_set = src_codon_set[1:i]
        codon_set_str = codon_set_to_str(codon_set)
        data = CodonGraphData(codon_set; plot_title = codon_set_str)
        # data_expanded = CodonGraphData(codon_set; plot_title = "$(string(codon_set)) after expanding")
        construct_graph_data!(data; debug = debug)
        # construct_graph_data!(data_expanded; debug = debug)
        push!(data_list, data)
    end
    show_multiple_codon_graphs(data_list; fig_title = "Line $line_count", debug = debug)
    return true
end

function grow_codon_results(src_codon_set::Vector{LongDNA{4}}; debug::Bool = false)
    for i in 1:length(src_codon_set)
        codon_set = src_codon_set[1:i]
        data = CodonGraphData(codon_set; plot_title = "$(string(codon_set)) before expanding")
        data_expanded = CodonGraphData(codon_set; plot_title = "$(string(codon_set)) after expanding")
        construct_graph_data!(data; debug = debug)
        construct_graph_data!(data_expanded; debug = debug)

        if is_strong_c3(data; debug = debug)
            data_list = [data, data_expanded]
            show_multiple_codon_graphs(data_list; debug = debug)
        else
            return false
        end
    end
    return true
end

function read_line(path, n::Int)
    open(path) do io
        for _ in 1:(n - 1)
            eof(io) && return nothing
            readline(io)
        end
        return eof(io) ? nothing : readline(io)
    end
end

function get_count_tier(codon::LongDNA{4}, size::Int)
    codon_count = get_codon_count_in_res("files/results/res_$(size).csv")
    count = codon_count[codon]
    count_of_count = get_count_of_count(codon_count)
    counts = sort(collect(keys(count_of_count)))

    min_count = counts[1]
    mid_count = counts[2]
    max_count = counts[3]

    if count == min_count
        return "min"
    elseif count == mid_count
        return "mid"
    elseif count == max_count
        return "max"
    else
        throw(
            ArgumentError(
                "Count $count does not match any tier (min: $min_count, mid: $mid_count, max: $max_count)",
            ),
        )
    end
end

function compare_analysis_files_linewise(orig_file_path::String, comp_file_path::String)
    open(orig_file_path, "r") do orig_file
        open(comp_file_path, "r") do comp_file
            line_number = 0
            for (orig_line, comp_line) in zip(eachline(orig_file), eachline(comp_file))
                line_number += 1
                if orig_line != comp_line
                    println("Difference found at line $line_number:")
                    println("Original: $orig_line")
                    println("Comparison: $comp_line")
                end
            end
        end
    end
end

# count which codons appear how many times in a results file
function get_codon_count_in_res(res_path::String)
    codon_counts = Dict{LongDNA{4}, Int}()

    # initialize counts to 0
    for codon in ALL_CODONS
        codon_counts[codon] = 0
    end

    open(res_path, "r") do f
        for line in eachline(f)
            codon_set = get_codon_set_from_res(line)
            for codon in codon_set
                codon_counts[codon] += 1
            end
        end
    end

    return codon_counts
end

function plot_codon_counts(res_path; order = ALL_CODONS, save_path = nothing)
    counts = get_codon_count_in_res(res_path)

    codons = order
    values = [counts[c] for c in codons]

    fig = Figure(size = (1800, 900))
    ax = Axis(fig[1, 1]; xlabel = "Codon", ylabel = "Anzahl", xticklabelrotation = pi / 3)

    barplot!(ax, 1:length(codons), values; color = :dodgerblue)
    ax.xticks = (1:length(codons), String.(codons))

    # counts als Labels oberhalb der Balken
    label_offset = maximum(values) == 0 ? 0.5 : maximum(values) * 0.01
    for (i, v) in enumerate(values)
        text!(
            ax,
            string(v);
            position = (i, v + label_offset),
            align = (:center, :bottom),
            fontsize = 10,
            color = :black,
        )
    end

    if save_path !== nothing
        save(save_path, fig)
    end

    return fig
end

function plot_codon_counts_sorted(res_path; order = ALL_CODONS, save_path = nothing)
    counts = get_codon_count_in_res(res_path)

    codons = order
    values = [counts[c] for c in codons]

    # nach Häufigkeit absteigend sortieren
    sorted_pairs = sort(collect(zip(codons, values)); by = last, rev = true)
    codons, values = first.(sorted_pairs), last.(sorted_pairs)

    fig = Figure(size = (1800, 900))
    ax = Axis(fig[1, 1]; xlabel = "Codon", ylabel = "Anzahl", xticklabelrotation = pi / 3)

    barplot!(ax, 1:length(codons), values; color = :dodgerblue)
    ax.xticks = (1:length(codons), String.(codons))

    # counts als Labels oberhalb der Balken
    label_offset = maximum(values) == 0 ? 0.5 : maximum(values) * 0.01
    for (i, v) in enumerate(values)
        text!(
            ax,
            string(v);
            position = (i, v + label_offset),
            align = (:center, :bottom),
            fontsize = 10,
            color = :black,
        )
    end

    if save_path !== nothing
        save(save_path, fig)
    end

    return fig
end

# function to get the amount of different values in codon_count dict
function get_count_of_count(codon_count::Dict{LongDNA{4}, Int})
    count_of_count = Dict{Int, Int}()
    for count in values(codon_count)
        if haskey(count_of_count, count)
            count_of_count[count] += 1
        else
            count_of_count[count] = 1
        end
    end
    return count_of_count
end

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
