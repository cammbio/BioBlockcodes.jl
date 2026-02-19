# validate CodonGraphData consistency (labels, indices, graph size)
function _validate_cgd(cgd::CodonGraphData)
    _validate_codon_set(cgd.codon_set)
    _validate_edges(cgd)
    _validate_vertices(cgd)
end


function _validate_ckp_file(path::AbstractString)
    # check if file exists
    !isfile(path) && throw(ArgumentError("checkpoint file not found.
    File path: \"$path\""))

    # check if file is empty
    (filesize(path) == 0) && throw(ArgumentError("checkpoint file is empty.
    File path: \"$path\""))

    # check if file has correct format (3 lines, each with key,value separated by comma)
    lines = readlines(path)
    length(lines) != 3 &&
        throw(ArgumentError("checkpoint file must contain exactly three lines with following format:
        comb_size,value1
        next_line,value2
        status,[finished|unfinished]
        With value1 1-20 and value2 1-n
        File path: \"$path\""))


    expected_keys = ("comb_size", "next_line", "status")
    vals_by_key = Dict{String, String}()
    # check if keys are valid
    for (line_count, line) in enumerate(lines)
        # check if line contains exactly one comma
        count(==(','), line) == 1 || throw(
            ArgumentError(
                "checkpoint file has invalid format at line $line_count. Expected exactly one comma. File path: \"$path\"",
            ),
        )
        # split line into key and value and validate strict formatting
        key, value = split(line, ","; limit = 2)

        # value must start directly after comma and no surrounding whitespaces are allowed
        (key == strip(key) && value == strip(value)) || throw(
            ArgumentError(
                "checkpoint file has invalid whitespace at line $line_count. Expected format \"key,value\" without spaces. File path: \"$path\"",
            ),
        )

        # check if key and value are not empty, key is allowed and not duplicated
        isempty(key) &&
            throw(ArgumentError("checkpoint file has empty key at line $line_count. File path: \"$path\""))
        isempty(value) &&
            throw(ArgumentError("checkpoint file has empty value at line $line_count. File path: \"$path\""))
        key == expected_keys[line_count] || throw(
            ArgumentError(
                "checkpoint file has invalid key \"$key\" at line $line_count. Expected \"$(expected_keys[line_count])\". File path: \"$path\"",
            ),
        )
        vals_by_key[key] = value
    end

    # check if values are valid
    try
        parse(Int, vals_by_key["comb_size"])
    catch
        throw(
            ArgumentError(
                "checkpoint file has non-integer value for \"comb_size\": \"$(vals_by_key["comb_size"])\". File path: \"$path\"",
            ),
        )
    end

    try
        parse(Int, vals_by_key["next_line"])
    catch
        throw(
            ArgumentError(
                "checkpoint file has non-integer value for \"next_line\": \"$(vals_by_key["next_line"])\". File path: \"$path\"",
            ),
        )
    end

    status = vals_by_key["status"]
    (status == "finished" || status == "unfinished") || throw(
        ArgumentError(
            "checkpoint file has invalid value for \"status\": \"$status\". Must be \"finished\" or \"unfinished\". File path: \"$path\"",
        ),
    )
end


function _validate_ckp(
    ckp::NamedTuple{(:comb_size, :next_line, :status)},
    ckp_path::String,
    comb_size::Int,
    max_lines::Int,
)
    ckp_comb_size = ckp.comb_size
    ckp_next_line = ckp.next_line
    ckp_status = ckp.status

    # check if checkpoint combination size matches the actual combination size
    (ckp_comb_size != comb_size) && throw(
        ArgumentError(
            "checkpoint file \"$ckp_path\" corrupted. Expected combination size: $comb_size, found: $ckp_comb_size.",
        ),
    )

    # check if checkpoint next line is valid
    (ckp_next_line < 0) &&
        throw(ArgumentError("checkpoint file \"$ckp_path\" corrupted. Invalid next line: $ckp_next_line."))
    (ckp_next_line > max_lines) && throw(
        ArgumentError(
            "checkpoint file \"$ckp_path\" corrupted. next line exceeds max lines: $ckp_next_line > $max_lines.",
        ),
    )

    # check if checkpoint status is valid
    (ckp_status != "finished" && ckp_status != "unfinished") && throw(
        ArgumentError(
            "checkpoint file \"$ckp_path\" corrupted. Invalid status: $ckp_status. Must be \"finished\" or \"unfinished\".",
        ),
    )

    # if checkpoint status is "finished" then next line must be 0
    (ckp_status == "finished" && ckp_next_line != 0) && throw(
        ArgumentError(
            "checkpoint file \"$ckp_path\" corrupted. Status is \"finished\" but next line \"$ckp_next_line\" is not \"0\", which indicates that status should be \"unfinished\".",
        ),
    )

    # if checkpoint status is "unfinished" then next line must not be 0
    (ckp_status == "unfinished" && ckp_next_line == 0) && throw(
        ArgumentError(
            "checkpoint file \"$ckp_path\" corrupted. Status is \"unfinished\" but next line \"$ckp_next_line\" is \"0\", which indicates that status should be \"finished\".",
        ),
    )
end


# validate codon_set contents
function _validate_codon(codon::LongDNA{4})
    # do not allow empty codon
    length(codon) == 0 && throw(ArgumentError("codon is empty."))

    # check if codon is valid
    codon in ALL_CODONS || throw(ArgumentError("codon \"$codon\" is not a valid codon."))
end


# validate codon_set contents
function _validate_codon_set(codon_set::Vector{LongDNA{4}})
    # do not allow empty codon sets
    length(codon_set) == 0 && throw(ArgumentError("\"codon_set\" is empty."))

    # do not allow duplicates
    length(codon_set) == length(Set(codon_set)) ||
        throw(ArgumentError("\"codon_set\" contains duplicate codons."))

    # check if all codons in codon_set are valid
    for codon in codon_set
        codon in ALL_CODONS || throw(ArgumentError("codon \"$codon\" in \"codon_set\" is not a valid codon."))
    end
end


function _validate_comb(comb::Vector{Int})
    # check if combination is empty
    length(comb) == 0 && throw(ArgumentError("\"comb\" is empty."))

    # check if combination contains valid indices
    max_index = length(ALL_CODONS)
    for index in comb
        (1 <= index <= max_index) ||
            throw(ArgumentError("\"comb\" contains invalid index: $index. Must be between 1 and $max_index."))
    end

    # check if combination contains duplicates
    length(comb) == length(Set(comb)) || throw(ArgumentError("\"comb\" contains duplicate indices."))

    # check if combination is sorted in ascending order
    issorted(comb) || throw(ArgumentError("\"comb\" must be sorted in ascending order."))
end


function _validate_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || throw(ArgumentError("directory does not exist: \"$dir\""))
end

# validate edge labels against edge count
function _validate_edges(cgd::CodonGraphData)
    # check if edge_labels is empty
    isempty(cgd.edge_labels) && throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" is empty."))
    # check if graph contains edges
    ne(cgd.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no edges."))
    # check if edge_labels length matches number of edges in graph
    length(cgd.edge_labels) == ne(cgd.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"edge_labels\" length does not match number of edges in graph.",
        ),
    )

    # check if duplicates exist in edge_labels
    length(cgd.edge_labels) == length(Set(cgd.edge_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" contains duplicate labels."))

    # check if every label in edge_labels corresponds to an actual edge in the graph
    for edge in cgd.edge_labels
        src_label, dst_label = edge
        # check if src_label exists in vert_idxs
        haskey(cgd.vert_idxs, src_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if dst_label exists in vert_idxs
        haskey(cgd.vert_idxs, dst_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if src_label exists in vert_labels
        src_label in cgd.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if dst_label exists in vert_labels
        dst_label in cgd.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if edge from src_label to dst_label exists in graph
        _has_edge_label(cgd, edge) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains edge from \"$src_label\" to \"$dst_label\" that does not exist in graph.",
            ),
        )
    end

    # check if every edge in graph has a corresponding label in edge_labels
    for edge in edges(cgd.graph)
        src_label = cgd.vert_labels[src(edge)]
        dst_label = cgd.vert_labels[dst(edge)]
        (src_label, dst_label) in cgd.edge_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: edge from \"$src_label\" to \"$dst_label\" does not have a corresponding label in \"edge_labels\".",
            ),
        )
    end
end


function _validate_graph(graph::SimpleDiGraph)
    # do not allow empty graphs
    ne(graph) == 0 && throw(ArgumentError("graph has no edges."))
    nv(graph) == 0 && throw(ArgumentError("graph has no vertices."))
end


function _validate_label(label::String)
    # do not allow empty labels
    isempty(label) && throw(ArgumentError("label cannot be empty."))
    # only allow labels of length 1 or 2
    (length(label) in (1, 2)) || throw(ArgumentError("label \"$label\" must be of length 1 or 2."))
    # only allow labels containing A, C, G and T
    for char in label
        char in ALLOWED_BASES_STR || throw(
            ArgumentError(
                "label \"$label\" contains invalid character \"$char\". Only A, C, G and T allowed.",
            ),
        )
    end
end

# validate vertices and index mappings
function _validate_vertices(cgd::CodonGraphData)
    # check if vert_labels is empty
    isempty(cgd.vert_labels) && throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" is empty."))
    # check if vert_idxs is empty
    isempty(cgd.vert_idxs) && throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" is empty."))
    # check if graph contains vertices
    nv(cgd.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no vertices."))
    # check if no duplicates exist in vert_labels
    length(cgd.vert_labels) == length(Set(cgd.vert_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" contains duplicate labels."))

    # check if vert_labels only contains valid labels (length 1 or 2, only A,C,G,T or combinations thereof)
    for label in cgd.vert_labels
        (length(label) in (1, 2)) || throw(
            ArgumentError("inconsistent CodonGraphData: \"vert_labels\" contains invalid label \"$label\"."),
        )
        for char in label
            char in ALLOWED_BASES_STR || throw(
                ArgumentError(
                    "inconsistent CodonGraphData: \"vert_labels\" contains invalid label \"$label\".",
                ),
            )
        end
    end

    # check if vert_labels length matches graph size
    length(cgd.vert_labels) == nv(cgd.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_labels\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs length matches graph size
    length(cgd.vert_idxs) == nv(cgd.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_idxs\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs keys match vert_labels and indices are valid
    for (label, idx) in cgd.vert_idxs
        (1 <= idx <= nv(cgd.graph)) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"vert_idxs\" contains invalid vertex index \"$idx\".",
            ),
        )
        cgd.vert_labels[idx] == label ||
            throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" and vert_labels do not match."))
    end
end
