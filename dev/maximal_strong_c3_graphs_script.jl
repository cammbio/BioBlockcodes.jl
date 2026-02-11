using Revise
using GCATCodes
using BioSequences
using Base.Threads
using CairoMakie

# stop_flag
const stop_flag = Base.Threads.Atomic{Bool}(false)


# grows a maximal strong c3 codon set from size 1-12 and displays graphs for each size
function get_cod_graphs(max_codon_set::Vector{LongDNA{4}}; debug::Bool = false)
    data_list = Vector{CodonGraphData}()
    for i in 1:length(max_codon_set)
        codon_set = max_codon_set[1:i]
        codon_set_str = codon_set_to_str(codon_set)
        data = CodonGraphData(codon_set; plot_title = codon_set_str)
        construct_graph_data!(data; debug = debug)
        push!(data_list, data)
    end
    return data_list
end


function run(max_res_path::String, stop_flag::Base.Threads.Atomic{Bool}; debug::Bool = false)
    stop_flag[] && return @info "Process was stopped before it started."

    open(max_res_path, "r") do io
        line_count = 1
        for line in eachline(io)
            if stop_flag[]
                debug && @debug "Process was stopped after processing $line_count lines."
                return false
            end
            codon_set = get_codon_set_from_res(line)
            data_list = get_cod_graphs(codon_set; debug = true)
            fig = show_multiple_codon_graphs(data_list; fig_title = "Line $line_count", debug = debug)

            save_path = joinpath(@__DIR__, "..", "files", "diagrams", "max", "line_$(string(line_count)).png")

            _save_graph(fig, save_path)
            fig = nothing
            # GC.gc()
            line_count += 1
        end
    end
    return true
end


function _save_graph(fig::Figure, save_path::String)
    save(save_path, fig)
end


println("Script loaded.")
