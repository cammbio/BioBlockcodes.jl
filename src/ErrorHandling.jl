# validate CodonGraphData consistency (labels, indices, graph size)
function _validate_cgd(data::CodonGraphData)
    _validate_codon_set(data.codon_set)
    _validate_edges(data)
    _validate_vertices(data)
end


function _validate_ckp_file(path::AbstractString)
    # 
    !isfile(path) && throw(ArgumentError("checkpoint file not found.
    File path: $path"))

    (filesize(path) == 0) && throw(ArgumentError("checkpoint file is empty.
    File path: $path"))

    countlines(path) != 3 &&
        throw(ArgumentError("checkpoint file must contain exactly three lines with following format:
        comb_size,value1
        curr_line,value2
        status,[finished|unfinished]
        With value1 1-20 and value2 1-n
        File path: $path"))
end


function _validate_checkpoint(ckp::Checkpoint, ckp_path::String, comb_size::Int, max_lines::Int)
    ckp_comb_size = ckp.comb_size
    ckp_curr_line = ckp.curr_line
    ckp_status = ckp.status

    # check if checkpoint combination size matches the actual combination size
    (ckp_comb_size != comb_size) && throw(
        ArgumentError(
            "Checkpoint file $ckp_path corrupted. Expected combination size: $comb_size, found: $ckp_comb_size.",
        ),
    )

    # check if checkpoint current line is valid
    (ckp_curr_line < 1) &&
        throw(ArgumentError("Checkpoint file corrupted. Invalid current line: $ckp_curr_line."))
end


# validate codon_set contents
function _validate_codon(codon::LongDNA{4})
    # do not allow empty codon
    length(codon) == 0 && throw(ArgumentError("\"codon\" is empty."))

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
    isdir(dir) || throw(ArgumentError("directory does not exist: $dir"))
end


# validate edge labels against edge count
function _validate_edges(data::CodonGraphData)
    # check if edge_labels is empty
    isempty(data.edge_labels) &&
        throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" is empty."))
    # check if graph contains edges
    ne(data.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no edges."))
    # check if edge_labels length matches number of edges in graph
    length(data.edge_labels) == ne(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"edge_labels\" length does not match number of edges in graph.",
        ),
    )

    # check if duplicates exist in edge_labels
    length(data.edge_labels) == length(Set(data.edge_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"edge_labels\" contains duplicate labels."))

    # check if every label in edge_labels corresponds to an actual edge in the graph
    for (src_label, dst_label) in data.edge_labels
        # check if src_label exists in vert_idxs
        haskey(data.vert_idxs, src_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if dst_label exists in vert_idxs
        haskey(data.vert_idxs, dst_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_idxs\".",
            ),
        )
        # check if src_label exists in vert_labels
        src_label in data.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$src_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if dst_label exists in vert_labels
        dst_label in data.vert_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains \"$dst_label\" that does not exist in \"vert_labels\".",
            ),
        )
        # check if edge from src_label to dst_label exists in graph
        _has_edge_label(data, src_label, dst_label) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"edge_labels\" contains edge from \"$src_label\" to \"$dst_label\" that does not exist in graph.",
            ),
        )
    end

    # check if every edge in graph has a corresponding label in edge_labels
    for edge in edges(data.graph)
        src_label = data.vert_labels[src(edge)]
        dst_label = data.vert_labels[dst(edge)]
        (src_label, dst_label) in data.edge_labels || throw(
            ArgumentError(
                "inconsistent CodonGraphData: edge from \"$src_label\" to \"$dst_label\" does not have a corresponding label in \"edge_labels\".",
            ),
        )
    end
end


# validate vertices and index mappings
function _validate_vertices(data::CodonGraphData)
    # check if vert_labels is empty
    isempty(data.vert_labels) &&
        throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" is empty."))
    # check if vert_idxs is empty
    isempty(data.vert_idxs) && throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" is empty."))
    # check if graph contains vertices
    nv(data.graph) == 0 && throw(ArgumentError("inconsistent CodonGraphData: graph has no vertices."))
    # check if no duplicates exist in vert_labels
    length(data.vert_labels) == length(Set(data.vert_labels)) ||
        throw(ArgumentError("inconsistent CodonGraphData: \"vert_labels\" contains duplicate labels."))

    # check if vert_labels only contains valid labels (length 1 or 2, only A,C,G,T or combinations thereof)
    for label in data.vert_labels
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
    length(data.vert_labels) == nv(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_labels\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs length matches graph size
    length(data.vert_idxs) == nv(data.graph) || throw(
        ArgumentError(
            "inconsistent CodonGraphData: \"vert_idxs\" length does not match number of vertices in graph.",
        ),
    )

    # check if vert_idxs keys match vert_labels and indices are valid
    for (label, idx) in data.vert_idxs
        (1 <= idx <= nv(data.graph)) || throw(
            ArgumentError(
                "inconsistent CodonGraphData: \"vert_idxs\" contains invalid vertex index \"$idx\".",
            ),
        )
        data.vert_labels[idx] == label ||
            throw(ArgumentError("inconsistent CodonGraphData: \"vert_idxs\" and vert_labels do not match."))
    end
end
