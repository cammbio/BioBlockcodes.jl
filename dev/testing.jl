using Revise
isdefined(Main, :BioBlockcodes) || using BioBlockcodes
using Base.Threads
using BioSequences
using CairoMakie
using Graphs
using GraphMakie
using JuliaFormatter
using Logging

# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
# global_logger(ConsoleLogger(Logging.Info)) # deactivate

const stop_flag = Base.Threads.Atomic{Bool}(false)
const amino_table = Dict(
    LongDNA{4}.("AAC") => "Asn",
    LongDNA{4}.("AAG") => "Lys",
    LongDNA{4}.("AAT") => "Asn",
    LongDNA{4}.("ACA") => "Thr",
    LongDNA{4}.("ACC") => "Thr",
    LongDNA{4}.("ACG") => "Thr",
    LongDNA{4}.("ACT") => "Thr",
    LongDNA{4}.("AGA") => "Arg",
    LongDNA{4}.("AGC") => "Ser",
    LongDNA{4}.("AGG") => "Arg",
    LongDNA{4}.("AGT") => "Ser",
    LongDNA{4}.("ATA") => "Ile",
    LongDNA{4}.("ATC") => "Ile",
    LongDNA{4}.("ATG") => "Met",
    LongDNA{4}.("ATT") => "Ile",
    LongDNA{4}.("CAA") => "Gln",
    LongDNA{4}.("CAC") => "His",
    LongDNA{4}.("CAG") => "Gln",
    LongDNA{4}.("CAT") => "His",
    LongDNA{4}.("CCA") => "Pro",
    LongDNA{4}.("CCG") => "Pro",
    LongDNA{4}.("CCT") => "Pro",
    LongDNA{4}.("CGA") => "Arg",
    LongDNA{4}.("CGC") => "Arg",
    LongDNA{4}.("CGG") => "Arg",
    LongDNA{4}.("CGT") => "Arg",
    LongDNA{4}.("CTA") => "Leu",
    LongDNA{4}.("CTC") => "Leu",
    LongDNA{4}.("CTG") => "Leu",
    LongDNA{4}.("CTT") => "Leu",
    LongDNA{4}.("GAA") => "Glu",
    LongDNA{4}.("GAC") => "Asp",
    LongDNA{4}.("GAG") => "Glu",
    LongDNA{4}.("GAT") => "Asp",
    LongDNA{4}.("GCA") => "Ala",
    LongDNA{4}.("GCC") => "Ala",
    LongDNA{4}.("GCG") => "Ala",
    LongDNA{4}.("GCT") => "Ala",
    LongDNA{4}.("GGA") => "Gly",
    LongDNA{4}.("GGC") => "Gly",
    LongDNA{4}.("GGT") => "Gly",
    LongDNA{4}.("GTA") => "Val",
    LongDNA{4}.("GTC") => "Val",
    LongDNA{4}.("GTG") => "Val",
    LongDNA{4}.("GTT") => "Val",
    LongDNA{4}.("TAA") => "Stop",
    LongDNA{4}.("TAC") => "Tyr",
    LongDNA{4}.("TAG") => "Stop",
    LongDNA{4}.("TAT") => "Tyr",
    LongDNA{4}.("TCA") => "Ser",
    LongDNA{4}.("TCC") => "Ser",
    LongDNA{4}.("TCG") => "Ser",
    LongDNA{4}.("TCT") => "Ser",
    LongDNA{4}.("TGA") => "Stop",
    LongDNA{4}.("TGC") => "Cys",
    LongDNA{4}.("TGG") => "Trp",
    LongDNA{4}.("TGT") => "Cys",
    LongDNA{4}.("TTA") => "Leu",
    LongDNA{4}.("TTC") => "Phe",
    LongDNA{4}.("TTG") => "Leu",
)


path = "files/results/sorted_res_12.csv"
path = "files/maximal_comma_free_strong_c3_codes/res_12.txt"
path = "files/maximal_self_complementary_strong_c3_codes/res_8.txt"
path = "files/maximal_self_complementary_comma_free_strong_c3_codes/res_8.txt"

amino_vector = Vector{Set{String}}()
for line in eachline(path)
    amino_set = Set{String}()
    codons = split(replace(line, "\"" => ""), ", ")
    codon_set = LongDNA{4}.(codons)
    cgd = CodonGraphData(codon_set)
    for codon in codon_set
        push!(amino_set, amino_table[codon])
    end
    push!(amino_vector, amino_set)
    println("Codon set: $codon_set, Amino acids: $amino_set, count of amino acids: $(length(amino_set))")
end

sorted_amino_vector = sort!(amino_vector; by = length)
for set in sorted_amino_vector
    println("Amino acid set: $set, count of amino acids: $(length(set))")
end

# get corresponding combination from codon set
function _get_comb_from_codon_set2(codon_set::Vector{LongDNA{4}})
    # get indices of codons in codon_set
    idxs = getindex.(Ref(BioBlockcodes.CODON_INDEX), codon_set)
    return idxs
end

function count_skipped(size::Int)
    len = 60
    path = "files/results/sorted_res_$(size-1).csv"
    total_combinations = binomial(len, size)
    generated_candidates = big(0)

    for line in eachline(path)
        codon_set = get_codon_set_from_line(line)
        comb = _get_comb_from_codon_set2(codon_set)
        last_idx = comb[end]
        generated_candidates += len - last_idx
    end

    skipped = total_combinations - generated_candidates
    return (; total_combinations, generated_candidates, skipped)
end

for i in 11:13
    sum_skipped = 0
    res = count_skipped(i)
    # println(
    #     "Size $i: Total combinations: $(res.total_combinations), Generated candidates: $(res.generated_candidates), Skipped: $(res.skipped)",
    # )
    # sum_skipped += res.skipped
    # println("Total skipped so far: $sum_skipped")
    println("Size: $i, Skipped: $(res.skipped)")
end

for i in 12:20
    println("Size $i: $(binomial(60, i))")
end
sum_skipped = BigInt(5166863418779)
for i in 14:20
    println("sum_skipped before size $i: $sum_skipped")
    println("Size $i: Skipped combinations: $(binomial(60, i))")
    count_skipped = binomial(60, i)
    sum_skipped += count_skipped
    println("sum_skipped after size $i: $sum_skipped")
end

sum_bin = big(0)
for i in 1:20
    bin = big(0)
    bin = binomial(60, i)
    println("Size $i: binomial(60, $i) = $bin")
    sum_bin += bin
    println("Sum of binomial coefficients up to size $i: $sum_bin")
end

a = big(7776048412324713)
b = big(1152921504606846975)
1 - (a / b)

binomial(60, 20)
codon_set1 = LongDNA{4}.(["ATG", "GCT"])
# codon_set2 = left_shift_codon_set(codon_set1, 1)
# codon_set3 = left_shift_codon_set(codon_set1, 2)

cgd = CodonGraphData(codon_set1, graph_title = "G(X)")
cgd2 = CodonGraphData(codon_set1, graph_title = "G'(X)")
_expand_graph!(cgd2)
# cgd3 = CodonGraphData(codon_set3, graph_title = "G(α₂(X))")
list = [cgd, cgd2]
plot_multiple_codon_graphs(list)


is_strong_c3(cgd)
is_comma_free(cgd)
is_circular(cgd)
is_c3(cgd)
plot_codon_graph(cgd)

# function to get all paths from a graph
function get_all_paths(data::CodonGraphData)
    g = data.graph
    labels = data.vert_labels
    verts = collect(vertices(g))
    paths = Vector{Vector{String}}()
    for s in verts, t in verts
        s == t && continue
        for p in all_simple_paths(g, s, t)
            push!(paths, labels[p])
        end
    end
    return paths
end

# function to get all paths from a graph that are longer than a given length
function get_paths_longer_than(data::CodonGraphData, min_len::Int)
    g = data.graph
    labels = data.all_vertex_labels
    verts = collect(vertices(g))
    paths = Vector{Vector{String}}()
    for s in verts, t in verts
        s == t && continue
        for p in all_simple_paths(g, s, t)
            length(p) > min_len && push!(paths, labels[p])
        end
    end
    return paths
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

function compare_files(orig_file_path::String, comp_file_path::String)
    open(orig_file_path, "r") do orig_file
        open(comp_file_path, "r") do comp_file
            line_number = 0
            for (orig_line, comp_line) in zip(eachline(orig_file), eachline(comp_file))
                line_number += 1
                if orig_line != comp_line
                    println("Difference found at line $line_number:")
                    println("Original:   $orig_line")
                    println("Comparison: $comp_line")
                    return false
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
            codon_set = _get_codon_set_from_line(line)
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

    return true
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

# function to get all subsets of a codon set
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


for i in 1:12
    println("Comparing files for combination size $i....................................................")
    res_path = "files/results/res_$(i).csv"
    ckp_path = "files/checkpoints/ckp_$(i).csv"
    sort_path = "files/results/sorted_res_$(i).csv"

    test_res_path = "files/tests/results/res_$(i).csv"
    test_ckp_path = "files/tests/checkpoints/ckp_$(i).csv"
    test_sort_path = "files/tests/results/sorted_res_$(i).csv"


    println("read 1: $(read(res_path) == read(test_res_path))")
    println("read 1: $(read(ckp_path) == read(test_ckp_path))")
    println("read 1: $(read(sort_path) == read(test_sort_path))")


    println("Countlines in res_$(i).csv:             $(countlines(res_path))")
    println("Countlines in test res_$(i).csv:        $(countlines(test_res_path))")
    println("Countlines in ckp_$(i).csv:             $(countlines(ckp_path))")
    println("Countlines in test ckp_$(i).csv:        $(countlines(test_ckp_path))")
    println("Countlines in sorted_res_$(i).csv:      $(countlines(sort_path))")
    println("Countlines in test sorted_res_$(i).csv: $(countlines(test_sort_path))")


    # println("Comparing res files:")
    # compare_files(res_path, test_res_path)
    println("Comparing ckp files:")
    compare_files(ckp_path, test_ckp_path)
    println("Comparing sorted res files:")
    compare_files(sort_path, test_sort_path)
end



for i in 12:12
    println("Checking i: $i")
    res_path = "D:/Bachelorarbeit_23_02_2026_11_35/BioBlockcodes/files/results/sorted_res_$i.csv"
    path_com_free = "files/maximal_comma_free_strong_c3_codes/strong_res_$i.txt"
    path_self_comp = "files/maximal_self_complementary_strong_c3_codes/strong_res_$i.txt"
    path_self_comp_com_free = "files/maximal_comma_free_self_complementary_strong_c3_codes/strong_res_$i.txt"

    com_free_cnt = 0
    self_comp_cnt = 0
    self_comp_com_free_cnt = 0

    open(path_com_free, "a") do comma_out
        open(path_self_comp, "a") do self_out
            open(path_self_comp_com_free, "a") do self_comma_out
                for line in eachline(res_path)
                    codon_set = get_codon_set_from_line(line)
                    codon_set_str = codon_set_to_str(codon_set)
                    cgd = CodonGraphData(codon_set)
                    _expand_graph!(cgd)
                    if is_comma_free(cgd)
                        com_free_cnt += 1
                        println(comma_out, codon_set_str)
                    end
                    if is_self_complementary(cgd)
                        self_comp_cnt += 1
                        println(self_out, codon_set_str)
                    end
                    if (is_comma_free(cgd)) && (is_self_complementary(cgd))
                        self_comp_com_free_cnt += 1
                        println(self_comma_out, codon_set_str)
                    end
                end
            end
        end
    end
    println("""for i $i:
            comma-free:           $com_free_cnt
            self_comp:            $self_comp_cnt
            self-comp-comma-free: $self_comp_com_free_cnt""")
end
