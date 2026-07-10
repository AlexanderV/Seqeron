---
type: concept
title: "Dot plot (word-match / k-tuple dot matrix)"
tags: [comparative-genomics, algorithm]
mcp_tools:
  - generate_dot_plot
sources:
  - docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
  - docs/algorithms/Comparative_Genomics/Dot_Plot_Generation.md
  - docs/Validation/reports/COMPGEN-DOTPLOT-001.md
source_commit: 37c54d6df345672ac216015eda5a1639544b4b01
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-dotplot-001-evidence
      evidence: "Test Unit ID: COMPGEN-DOTPLOT-001 ... Algorithm: Dot Plot Generation (word-match / k-tuple dot matrix)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:average-nucleotide-identity
      source: compgen-dotplot-001-evidence
      evidence: "Dot plot (COMPGEN-DOTPLOT-001) and ANI (COMPGEN-ANI-001) are sibling Comparative-genomics pairwise-comparison units — the dot plot is the visual word-match matrix, ANI the averaged fragment-identity metric of the same two-sequence comparison"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:synteny-and-rearrangement-detection
      source: compgen-dotplot-001-evidence
      evidence: "Wikipedia dot plot: dots combine into diagonal lines (similarity), indels disrupt the diagonal, repeats add extra diagonals — the classic visual read-out of synteny and rearrangement between two sequences"
      confidence: medium
      status: current
---

# Dot plot (word-match / k-tuple dot matrix)

A **dot plot** is the classic two-dimensional visual comparison of two sequences: a matrix
whose axes are the two sequences, with a **dot at position (i, j)** wherever a **word (k-tuple)
of length `wordSize` starting at position i in sequence 1 exactly matches** the word starting
at position j in sequence 2. Regions of similarity appear as **diagonal runs of dots**;
indels shift/break a diagonal, repeats produce extra parallel diagonals, and a
self-comparison always shows the full **main diagonal**. It is a Comparative-genomics
family (`COMPGEN-*`) unit, a sibling of [[average-nucleotide-identity]],
[[synteny-and-rearrangement-detection]],
[[conserved-gene-clusters-common-intervals]], and the end-to-end pipeline
[[genome-comparison-core-dispensable]]: where ANI reduces the comparison to one *number* and
synteny to *ordered blocks*, the dot plot keeps the **whole match relation** as a visual
matrix. Validated under test unit **COMPGEN-DOTPLOT-001**; the pre-implementation evidence is
[[compgen-dotplot-001-evidence]] and the independent two-stage verdict is
[[compgen-dotplot-001-report]] (Stage A PASS / Stage B PASS-WITH-NOTES / End state CLEAN — no code
defect; one weak self-diagonal assertion strengthened to the exact main-diagonal set and the
default `wordSize=10` path newly covered, both in-session). [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

## The word-match algorithm (Gibbs & McIntyre 1970 → EMBOSS `dottup`)

The dot-matrix ("diagram") method was introduced by **Gibbs & McIntyre 1970**. Seqeron
implements the **word-match** variant — the exact-match form embodied by EMBOSS `dottup`
(*"looks for places where words (tuples) of a specified length have an exact match in both
sequences and draws a diagonal line over the position of these words"*), **not** the
substitution-matrix scoring form (`dotmatcher`).

1. For every word-start position `x` in sequence 1 (the x-axis) and `y` in sequence 2 (the
   y-axis), the length-`wordSize` words are compared.
2. A **dot (x, y)** is emitted exactly when those words are **equal** (case-insensitive; both
   sequences are upper-cased). At `wordSize = 1` this reduces to *"place a dot where the two
   characters are equal"* — i.e. a k-mer match with k = 1.
3. All (overlapping) occurrences are reported — a word repeated in sequence 2 yields a dot for
   each occurrence.

The exact-match engine finds occurrences via a suffix tree, so a word appearing multiple times
contributes multiple dots on its row — the same generalized-suffix-tree machinery behind
[[longest-common-substring]] (which instead returns the single deepest two-string node).

### Word size — the sensitivity/noise trade-off

The single most important parameter is **`wordSize`** (EMBOSS default **10**):

- **Longer words** → less random noise, faster, but **less sensitive** (miss short/fragmentary
  similarities).
- **Shorter words** → more sensitive to short similar regions, but **more spurious off-diagonal
  dots** (chance matches) and slower.

This is a documented trade-off, **not a defect**: the probability of matching several residues
in a row by chance is much lower than a single-residue match, which is precisely why tuples
suppress noise relative to `wordSize = 1`.

## Parameters and invariants

| Parameter | Meaning | Default |
|-----------|---------|---------|
| `wordSize` | length of the exact-match word (k-tuple) | 10 (EMBOSS `dottup`) |
| `stepSize` | sample every `stepSize`-th x-axis word-start position | 1 |

- **Coordinate orientation**: sequence 1 → x, sequence 2 → y. (An *assumption* — the choice of
  which input is the x-axis is a presentation convention; transposing would mirror the plot but
  not change the match set as a relation.)
- **Case-insensitive**: `a` matches `A` (both upper-cased); standard for nucleotide/protein dot
  plots.
- **Self-comparison** (sequence1 = sequence2) always yields the **full main diagonal**
  (`(i, i)` for all i at `wordSize = 1`), plus any off-diagonal chance matches.
- Number of dots ≤ (word-start positions) × (occurrences) — the `O(n × m)` complexity bound.

## Documented oracles

- **k = 1 worked example** (Huttley): sequence1 `AGCGT` (x), sequence2 `AT` (y), `wordSize = 1`
  → exactly `{(0,0)` (A = A)`, (4,1)` (T = T)`}`.
- **Exact word match, `wordSize = 4`** (EMBOSS `dottup`, all overlapping occurrences):
  `ACGTACGT` vs itself → every length-4 word start x = 0..4 →
  `{(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)}` (word `ACGT` at x = 0 and x = 4 each match at
  y = 0 and y = 4).
- **Self-comparison main diagonal**: `ACGT` vs itself, `wordSize = 1` → dots include
  `(0,0),(1,1),(2,2),(3,3)` (the full main diagonal) plus any off-diagonal chance matches.
- **`stepSize = 2`** samples every other x → only even x coordinates appear.

## Corner cases and failure modes

- **Word longer than a sequence** → no word of length `wordSize` can be formed → **no dots**.
- **Null / empty** input → no dots.
- **Completely dissimilar sequences** (disjoint alphabets) → **no dots** (a dot is drawn only on
  a match).
- **Random noise with small words** → spurious off-diagonal dots; this is the documented
  sensitivity/noise trade-off, not a bug.
- **Non-positive `wordSize` / `stepSize`** → `ArgumentOutOfRangeException` (a non-positive window
  is undefined for `dottup`).

## Reference tools

The method traces to **Gibbs & McIntyre 1970** (*Eur. J. Biochem.* 16(1):1–11, the originating
dot-matrix "diagram" method), **EMBOSS `dottup`** (Rice, Longden & Bleasby 2000 — the word-match
reference implementation, default `wordsize` 10), **Wikipedia — Dot plot (bioinformatics)** (the
dot-placement rule, diagonals-as-similarity, main-diagonal, and tuple noise-reduction), and the
**Huttley TIB Dotplot** course material (the k = 1 worked example). No deviations from the
sources are recorded; the only assumptions are the x = sequence1 orientation convention and
case-insensitive comparison, neither correctness-affecting.
