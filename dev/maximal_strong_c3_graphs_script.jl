using Revise
using BioSequences: LongDNA
using CairoMakie
using BioBlockcodes

# stop_flag
const stop_flag = Base.Threads.Atomic{Bool}(false)


function run(max_res_path::String, stop_flag::Base.Threads.Atomic{Bool})
    stop_flag[] && return false

    open(max_res_path, "r") do io
        line_count = 1
        for line in eachline(io)
            if stop_flag[]
                println("Process was stopped after processing $line_count lines.")
                return false
            end
            codon_set = get_codon_set_from_line(line)
            cgd_list = _get_cod_graphs(codon_set)
            fig = show_multiple_codon_graphs(cgd_list; fig_title = "Line $line_count")

            save_path = joinpath(@__DIR__, "..", "files", "diagrams", "max", "line_$(string(line_count)).png")

            _save_graph(fig, save_path)
            fig = nothing
            line_count += 1
        end
    end
    return true
end


# grows a maximal strong c3 codon set from size 1-12 and displays graphs for each size
function _get_cod_graphs(max_codon_set::Vector{LongDNA{4}})
    cgd_list = Vector{CodonGraphData}()
    for i in 1:length(max_codon_set)
        codon_set = max_codon_set[1:i]
        codon_set_str = codon_set_to_str(codon_set)
        cgd = CodonGraphData(codon_set; graph_title = codon_set_str)
        push!(cgd_list, cgd)
    end
    return cgd_list
end


function _save_graph(fig::Figure, save_path::String)
    save(save_path, fig)
end


println("Script loaded.")
