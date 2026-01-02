# ---------------------------------------------- VARIABLES ----------------------------------------------

# ---------------------------------------------- CONSTANTS ----------------------------------------------

# ---------------------------------------------- FUNCTIONS ----------------------------------------------
"""
    display_cycles(data::CodonGraphData; show_debug::Bool = false)

Find and display all cycles in `data.graph`.

# Arguments

  - `data::CodonGraphData`: Graph data to analyze.

# Keyword Arguments

  - `show_debug::Bool`: Whether to emit debug logs.

# Returns

  - `Nothing`: Prints cycle info when `show_debug` is true.

# Throws

  - None.

# Example

```julia
data = CodonGraphData(LongDNA{4}.(["CGT", "GTA", "ACT", "AAT"]))
construct_graph_data!(data)
display_cycles(data; show_debug = true)
```
"""
function display_cycles(data::CodonGraphData; show_debug::Bool = false)
    cycles = simplecycles(data.graph)
    show_debug && @debug """Found $(length(cycles)) cycles in graph.
    Cycles: $cycles"""

    # check for duplicates
    keys = map(Tuple, cycles)
    unique_keys = unique(keys)
    show_debug && @debug """Unique cycles: $unique_keys
    Total cycles: $keys"""
    if length(keys) != length(unique_keys)
        show_debug && @debug "Duplicate cycles found!"
        return false
    end

    # iterate all cycles and print them
    for cycle in cycles # iterate all cycles
        # get all vertice labels 
        # join every vertice label in the cycle with " -> " and print it
        show_debug && join((data.all_vertex_labels[i] for i in (cycle..., first(cycle))), " -> ")
        show_debug && @debug "Cycle length: $(length(cycle))"
    end
    show_debug && @debug "Amount of cycles found: $(length(cycles))"
end
