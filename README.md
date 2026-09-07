# `BioBlockcodes.jl`

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cammbio.github.io/BioBlockcodes.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cammbio.github.io/BioBlockcodes.jl/dev/)

`BioBlockcodes.jl` is a package to analyze circular codes [1].

A set of codons (words of length 3 over the DNA alphabet) is a *circular code* if every sequence written on a circle can be split into codons in exactly one way, no matter where the reading starts. Circular codes are studied in theoretical biology because they let the reading frame be recovered locally, without a start signal.

`BioBlockcodes.jl` decides these properties with a graph criterion: it builds the *codon graph* of a codon set and checks it for cycles. On top of this it provides predicates for a codon set

  - `is_circular` – the codon graph is acyclic;
  - `is_c3` – the code and both of its reading-frame shifts are circular;
  - `is_strong_c3` – C3 and the expanded codon graph has no cycle longer than 2;

as well as further checks such as `is_comma_free` and `is_self_complementary`. Each predicate accepts either a raw `Vector{LongDNA}` codon set or a prebuilt `CodonGraphData` object, so the graph can be reused across several tests. An optional package extension, `BioBlockcodesGraphMakie`, draws the codon graphs once a [Makie](https://docs.makie.org) backend and `GraphMakie` are loaded.

See the [tutorial](https://cammbio.github.io/BioBlockcodes.jl/dev/) for a step-by-step walkthrough of building codon sets, running the predicates, and plotting the resulting codon graphs.

`BioBlockcodes.jl` is the Julia version of the GCAT (**G**enetic **C**ode **A**nalysis **T**oolkit) written in Java [2].


## References

1 Fimmel, Elena, und Lutz Strüngmann. 2018. „Mathematical Fundamentals for the Noise Immunity of the Genetic Code“. Biosystems 164 (Februar): 186–98. [https://doi.org/10.1016/j.biosystems.2017.09.007](https://doi.org/10.1016/j.biosystems.2017.09.007).

2 Kraljić, K., L. Strüngmann, E. Fimmel, und M. Gumbel. 2018. „Genetic Code Analysis Toolkit: A Novel Tool to Explore the Coding Properties of the Genetic Code and DNA Sequences“. SoftwareX 7 (Januar): 12–14. [https://doi.org/10.1016/j.softx.2017.10.008](https://doi.org/10.1016/j.softx.2017.10.008).