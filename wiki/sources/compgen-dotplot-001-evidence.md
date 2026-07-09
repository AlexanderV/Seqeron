---
type: source
title: "Evidence: COMPGEN-DOTPLOT-001 (Dot plot — word-match / k-tuple)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
source_commit: e932a6407e97bc684626ab5e2012dccf91053200
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-DOTPLOT-001

The validation-evidence artifact for test unit **COMPGEN-DOTPLOT-001** — **Dot Plot
Generation** (word-match / k-tuple dot matrix): emit a dot at `(x, y)` wherever a
length-`wordSize` word starting at x in sequence 1 exactly matches the word starting at y in
sequence 2. This is a **Comparative-genomics** family Evidence file and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm, its parameters, invariants, worked oracles, and corner cases are summarized in
[[dot-plot-word-match]]. Its sibling COMPGEN units are [[average-nucleotide-identity]],
[[synteny-and-rearrangement-detection]], and
[[conserved-gene-clusters-common-intervals]]. See [[test-unit-registry]] for how units are
tracked.

## What this file records

- **Online sources:**
  - **Gibbs & McIntyre 1970** (*Eur. J. Biochem.* 16(1):1–11, DOI 10.1111/j.1432-1033.1970.tb01046.x,
    authority rank 1) — the originating dot-matrix ("diagram") method for comparing two
    sequences. The full text is paywalled (HTTP 403), so the method description is taken from
    the secondary sources that cite it; only the citation/DOI is attributed to this primary.
  - **Wikipedia — Dot plot (bioinformatics)** (rank 4, citing Gibbs & McIntyre 1970) — the
    dot-placement rule (a dot at row i / column j means residues at i and j match), diagonals =
    similarity ("identical proteins ... a diagonal line in the center"), main diagonal =
    self-alignment, and noise reduction via tuples (shading runs of k residues because the
    chance of k-in-a-row matching is much lower).
  - **EMBOSS `dottup` manual + manpage** (rank 3, reference implementation) — word-match dot
    plot: "looks for places where words (tuples) of a specified length have an exact match in
    both sequences" (exact matching, not `dotmatcher` scoring); the longer-word/less-noise vs
    shorter-word/more-sensitive trade-off; default `wordsize` **10** (the source of the impl's
    `wordSize` default).
  - **Huttley — TIB Dotplot** (rank 4, course material) — the k = 1 rule (`if X[i] == Y[j] then
    match`, "a k-mer matching algorithm where k = 1") and the `AGCGT` vs `AT` worked example
    → dots at (0,0) and (4,1).
- **Algorithm behaviour (from the artifact):** for each word-start x in seq1 and y in seq2,
  emit dot (x, y) when the length-`wordSize` words are equal (case-insensitive, both
  upper-cased); all overlapping occurrences reported (suffix-tree lookup); `stepSize` samples
  every stepSize-th x. x = sequence1, y = sequence2.
- **Datasets (documented oracles):**
  - *Huttley k = 1*: `AGCGT` (x) vs `AT` (y), wordSize 1 → exactly {(0,0), (4,1)}.
  - *Exact word match, wordSize 4*: `ACGTACGT` vs itself → every word start x = 0..4 →
    {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} (all overlapping occurrences).
  - *Self-comparison main diagonal*: `ACGT` vs itself, wordSize 1 → dots include the full main
    diagonal (0,0),(1,1),(2,2),(3,3) plus any off-diagonal chance matches.

## Deviations and assumptions

**Deviations: none.** Two **assumptions**, both explicitly non-correctness-affecting:

1. **Coordinate orientation** (x = sequence1, y = sequence2) — the sources fix the dot at
   (position-in-A, position-in-B) but which input is the x-axis is a presentation convention;
   transposing mirrors the plot without changing the match set as a relation.
2. **Case-insensitive comparison** — `dottup`/Gibbs do not mandate case folding; the impl
   upper-cases both sequences so `a` matches `A`, standard for nucleotide/protein dot plots
   where case is not biologically meaningful (documented in the algorithm doc).

Recommended coverage (MUST): k = 1 example {(0,0),(4,1)}; wordSize 4 self-match full overlapping
set; self-comparison main diagonal; disjoint-alphabet → no dots; word longer than sequence /
null / empty → no dots; non-positive wordSize/stepSize → `ArgumentOutOfRangeException`. SHOULD:
stepSize 2 samples even x only; case-insensitivity. No contradictions among sources — Gibbs &
McIntyre (via secondaries), Wikipedia, EMBOSS `dottup`, and Huttley agree on the exact-word
match rule, diagonals-as-similarity, and the wordSize noise/sensitivity trade-off.
