```@meta
CurrentModule = BioBlockcodes
```
Documentation for [BioBlockcodes](https://github.com/cammbio/BioBlockcodes.jl).

# Tutorial

First of all, load the required packages:

```jldoctest is_circular_tutorial
julia> using BioBlockcodes, BioSequences
```

`BioSequences` provides the `dna"..."` string macro that produces `LongDNA{4}` values which are required later.

## Testing a codon set for circularity with `is_circular`

This tutorial walks you through `is_circular(codons::Vector{LongDNA{4}})`, the
convenience method that answers a single question: *do the given codons form a
circular code?*

### Background

A set of codons (words of length 3 over the DNA alphabet) is a **circular code**
if every sequence written on a circle can be decoded into codons in only one way,
regardless of where you start reading. Circular codes are of interest in
theoretical biology because they allow the reading frame to be retrieved locally,
without a start signal.

`BioBlockcodes` decides circularity with a graph criterion: it builds the
*codon graph* of the set (vertices are the mononucleotides and dinucleotides
occurring in the codons, and each codon `N1N2N3` contributes the edges
`N1 -> N2N3` and `N1N2 -> N3`) and checks that this directed graph contains no
cycle. A codon set is circular **if and only if** its codon graph is acyclic.

### The two methods of `is_circular`

`is_circular` is defined for two argument types:

  - [`is_circular(cgd::CodonGraphData)`](@ref) works on a codon graph you have
    already constructed.
  - `is_circular(codons::Vector{LongDNA{4}})` takes the raw codon set, builds the
    [`CodonGraphData`](@ref) for you, and forwards to the method above.

Use the vector method when you just want a yes/no answer and do not need the
graph object for anything else.


### Build a codon set

A codon set is simply a `Vector{LongDNA{4}}` in which every entry has length 3:

```jldoctest is_circular_tutorial
julia> codon_set = [dna"AAC", dna"GTT"]
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 AAC
 GTT
```

You can also build the vector from strings, which is handy for longer sets:

```jldoctest is_circular_tutorial
julia> LongDNA{4}.(["AAC", "GTT"])
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 AAC
 GTT
```

### Call `is_circular`

```jldoctest is_circular_tutorial
julia> is_circular(codon_set)
true
```

The codon graph of `{AAC, GTT}` has no directed cycle, so the set is a circular
code.

### A non-circular example

The set `{AGA, AAG, GAA}` consists of the three cyclic permutations of the same
pattern. Its codon graph contains a cycle, so it is **not** a circular code:

```jldoctest is_circular_tutorial
julia> is_circular([dna"AGA", dna"AAG", dna"GAA"])
false
```

### Relationship to the `CodonGraphData` method

The following two calls are equivalent; the vector method just spares you the
explicit construction step:

```jldoctest is_circular_tutorial
julia> is_circular([dna"AAC", dna"GTT"])
true

julia> cgd = CodonGraphData([dna"AAC", dna"GTT"]);

julia> is_circular(cgd)
true
```

Reach for the `CodonGraphData` method when you also need to run other checks on
the same set (for example [`is_c3`](@ref), [`is_comma_free`](@ref) or
[`is_self_complementary`](@ref)), so the graph is built only once.

### Notes and errors

  - Every codon in `codons` must be a valid length-3 DNA sequence. An invalid
    codon set (wrong length, empty, duplicates, …) raises an `ArgumentError`
    while the `CodonGraphData` is being constructed.
  - The order of the codons in the vector does not affect the result.
  - `is_circular` tests circularity only. Stronger properties such as
    comma-freeness or the C3 property have their own predicates.

# API 

```@index
```

```@autodocs
Modules = [BioBlockcodes]
```
