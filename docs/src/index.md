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

## Testing a codon set for the C3 property with `is_c3`

This section covers `is_c3(codons::Vector{LongDNA{4}})`, which decides whether a
codon set is a *C3 code*. The setup and codon-set construction from the
`is_circular` section above carry over unchanged; only the property being tested
is new.

### Background

A codon set is **C3** if the set itself *and* both of its reading-frame shifts
are circular. The two frame shifts are obtained by rotating every codon
`N1N2N3` cyclically by one and by two positions; `is_c3` builds the codon graph
for each of the three sets and checks that none of them contains a cycle.

C3 sits strictly between the two other predicates in this tutorial: every C3 code
is [`is_circular`](@ref), and every [`is_strong_c3`](@ref) code is C3, but
neither implication reverses.

### The two methods of `is_c3`

As with `is_circular`, the predicate is defined for both argument types:

  - [`is_c3(cgd::CodonGraphData)`](@ref) works on an existing codon graph.
  - `is_c3(codons::Vector{LongDNA{4}})` builds the [`CodonGraphData`](@ref)
    internally and forwards to the method above. Note that the two frame-shifted
    sets get their own graphs regardless of which method you call.

### Checking a codon set

```jldoctest is_circular_tutorial
julia> is_c3([dna"AAC", dna"ATG"])
true
```

`{AAC, ATG}` is circular and stays circular under both frame shifts, so it is a
C3 code.

### Circular but not C3

Circularity of the set alone is not enough. `{AAC, CCA}` is circular, but one of
its frame shifts is not, so it fails the C3 test:

```jldoctest is_circular_tutorial
julia> is_circular([dna"AAC", dna"CCA"])
true

julia> is_c3([dna"AAC", dna"CCA"])
false
```

### Relationship to the `CodonGraphData` method

The two calls below are equivalent; pass a `CodonGraphData` when you want to run
several checks on the same set without rebuilding its graph:

```jldoctest is_circular_tutorial
julia> is_c3([dna"AAC", dna"ATG"])
true

julia> cgd = CodonGraphData([dna"AAC", dna"ATG"]);

julia> is_c3(cgd)
true
```

### Notes and errors

  - An invalid codon set raises an `ArgumentError` while the `CodonGraphData` is
    being constructed, exactly as for `is_circular`.
  - The order of the codons does not affect the result.
  - A `false` result does not say which frame is to blame; call
    [`is_circular`](@ref) on the set to see whether circularity already fails or
    only one of the shifts does.

## Testing a codon set for the strong C3 property with `is_strong_c3`

This section covers `is_strong_c3(codons::Vector{LongDNA{4}})`, which decides
whether a codon set is a *strong comma-free C3 code*. The setup and codon-set
construction from the `is_circular` section above carry over unchanged; only the
property being tested is new.

### Background

The relevant properties form a hierarchy, each one strictly stronger than the
one before:

  1. **circular** – the codon graph is acyclic (see the previous section).
  2. **C3** – the code *and* both of its frame shifts (obtained by cyclically
     shifting every codon one and two positions to the left) are circular. This
     is what [`is_c3`](@ref) checks.
  3. **strong C3** – the code is C3 *and* its *expanded* codon graph contains no
     cycle longer than length 2. The expanded graph adds, for every codon
     `N1N2N3`, the vertices `N2` and `N3N1` together with the two edges
     `N2 -> N3N1` and `N3N1 -> N2`.

So `is_strong_c3(codons)` returns `true` only for sets that already pass
`is_circular` and `is_c3`.

### The two methods of `is_strong_c3`

Like `is_circular`, the predicate is defined for both argument types:

  - [`is_strong_c3(cgd::CodonGraphData)`](@ref) works on an existing codon graph.
  - `is_strong_c3(codons::Vector{LongDNA{4}})` builds the [`CodonGraphData`](@ref)
    internally and forwards to the method above.

### Checking a codon set

```jldoctest is_circular_tutorial
julia> is_strong_c3([dna"GGA", dna"TAA"])
true
```

The set `{GGA, TAA}` is circular, stays circular under both frame shifts, and its
expanded graph has only length-2 cycles, so it is strong C3.

### C3 but not strong C3

Passing `is_c3` is necessary but not sufficient. The set `{AAC, ACC}` is C3, yet
its expanded graph contains a longer cycle, so it fails the strong C3 test:

```jldoctest is_circular_tutorial
julia> is_c3([dna"AAC", dna"ACC"])
true

julia> is_strong_c3([dna"AAC", dna"ACC"])
false
```

### A codon set that fails outright

`{AGA, AAG, CTG, TGA, TTC}` is not even circular, so every stronger property
fails as well:

```jldoctest is_circular_tutorial
julia> is_strong_c3([dna"AGA", dna"AAG", dna"CTG", dna"TGA", dna"TTC"])
false
```

### Relationship to the `CodonGraphData` method

As with `is_circular`, the two calls below are equivalent; pass a
`CodonGraphData` when you want to run several checks on the same set without
rebuilding the graph each time:

```jldoctest is_circular_tutorial
julia> is_strong_c3([dna"GGA", dna"TAA"])
true

julia> cgd = CodonGraphData([dna"GGA", dna"TAA"]);

julia> is_strong_c3(cgd)
true
```

### Notes and errors

  - An invalid codon set raises an `ArgumentError` while the `CodonGraphData` is
    being constructed, exactly as for `is_circular`.
  - The order of the codons does not affect the result.
  - `is_strong_c3` returning `false` tells you the set is not strong C3, but not
    which weaker property it still satisfies; call [`is_c3`](@ref) or
    [`is_circular`](@ref) to locate the set in the hierarchy.

# API 

```@index
```

```@autodocs
Modules = [BioBlockcodes]
```
