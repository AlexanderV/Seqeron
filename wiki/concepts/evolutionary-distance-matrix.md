---
type: concept
title: "Evolutionary distance matrix (p-distance / Jukes–Cantor / Kimura-2-parameter / Hamming — the substrate distance-based tree building consumes)"
tags: [phylogenetics, algorithm]
sources:
  - docs/Evidence/PHYLO-DIST-001-Evidence.md
source_commit: 3a53115ec5fbdbc54448d69550c3b961c40a320a
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: phylo-dist-001-evidence
      evidence: "Test Unit ID: PHYLO-DIST-001 ... Title: Phylogenetic Distance Matrix Calculation"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:phylogenetic-bootstrap-support
      source: phylo-dist-001-evidence
      evidence: "PHYLO-BOOT-001 rebuilds each replicate tree by the distance-matrix → UPGMA/NJ machinery; its oracles use UPGMA + JukesCantor distances (this unit) — bootstrap wraps this distance step"
      confidence: high
      status: current
---

# Evolutionary distance matrix (p-distance / Jukes–Cantor / Kimura-2-parameter / Hamming)

The **third phylogenetics-family (`PHYLO-*`) unit**, PHYLO-DIST-001 — the **pairwise evolutionary
distance** primitive: given a set of **aligned** sequences, produce the symmetric *n×n* matrix of
pairwise distances that **distance-based tree construction (Neighbor-Joining / UPGMA) consumes**. This
is the substrate underneath the rest of the family: [[phylogenetic-bootstrap-support]] (PHYLO-BOOT-001)
rebuilds every replicate tree by the same **distance-matrix → UPGMA/NJ** machinery (its oracles use
`UPGMA` + `JukesCantor` distances), and [[tree-comparison-metrics]] (PHYLO-COMP-001) compares the
`PhyloNode` trees that come *out* of that builder. Validated under test unit **PHYLO-DIST-001**; the
literature-traced record is [[phylo-dist-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade correctness
reference ([[scientific-rigor|research-grade]]), not for clinical use.

The public surface is two methods: `CalculatePairwiseDistance(s1, s2, method)` for a single pair and
`CalculateDistanceMatrix(seqs, method)` for the whole *n×n* matrix, over four `method` values —
**Hamming**, **PDistance**, **JukesCantor**, **Kimura2Parameter**.

## The four distance methods

All four operate **column-wise over aligned positions** and count/ratio *substitutions* against the
number of **comparable sites** (positions where both sequences have a real nucleotide — gaps and
ambiguous IUPAC bases are skipped; see below).

- **Hamming** — the **raw integer count** of differing positions. Not a proportion; the only method
  that does not divide by comparable sites.
- **p-distance** (Hamming *proportion*) — `p = differences / comparableSites` ∈ [0,1]. The observed
  fraction of mismatched sites, with **no correction** for unseen (multiple/back) substitutions.
- **Jukes–Cantor 1969 (JC69)** — `d = -3/4 · ln(1 − 4p/3)`, correcting p for multiple hits at the same
  site under the **JC model** (equal base frequencies π = 0.25 each, equal rates between all
  nucleotides). Always `d ≥ p`.
- **Kimura 2-parameter (K80/K2P)** — `d = -1/2 · ln((1 − 2S − V) · √(1 − 2V))`, where **S** = proportion
  of **transition** differences (purine↔purine A↔G, pyrimidine↔pyrimidine C↔T) and **V** = proportion of
  **transversion** differences. K2P separates the transition/transversion rate; also `d ≥ p`.

## Matrix properties and invariants (§4)

1. **Symmetry** — `M[i,j] = M[j,i]` for all *i,j* (time-reversible models; `d(s1,s2)=d(s2,s1)`).
2. **Zero diagonal** — `M[i,i] = 0` (a sequence's distance to itself).
3. **Non-negative** — `M[i,j] ≥ 0` for the corrected distances.
4. **Dimensions** — an *n×n* matrix for *n* input sequences.
5. **Correction ordering** — `JC69 ≥ p-distance` and `K2P ≥ p-distance`: the model corrections only
   ever *inflate* the raw proportion (they estimate the hidden substitutions p cannot see).

The **triangle inequality** `d(A,C) ≤ d(A,B) + d(B,C)` is *expected to hold in most cases* but is **not
guaranteed** by the corrected (non-metric) distances — the source flags this explicitly rather than
asserting a metric.

## Gap and ambiguity handling (pairwise deletion)

- **Gaps (`-`)** are **excluded** from the comparison — a position with a gap in *either* sequence is
  skipped and is **not** counted in `comparableSites`. (Gaps-as-mismatch is a documented alternative in
  the literature; this unit ignores them.)
- **Ambiguous IUPAC bases** (`N`, `R`, `Y`, …) are **skipped like gaps** — only A, C, G, T are
  compared (pairwise deletion).
- **Case-insensitive** — upper/lowercase are treated the same (`ToUpperInvariant`), standard
  bioinformatics practice.

## Saturation and corner cases (§3, §2)

The correction formulas have a **domain limit** — beyond it the log argument goes non-positive and the
distance is **undefined / +∞** (biological *saturation*: too many substitutions to estimate):

- **JC69 saturates at `p ≥ 3/4`** — `1 − 4p/3 ≤ 0` → `+∞`. (Two maximally divergent sequences over a
  4-letter alphabet share ~25% of sites by chance.)
- **K2P saturates at `V ≥ 1/2`** — `1 − 2V ≤ 0` → `+∞`.

Other documented corner cases:

| Case | Behavior |
|------|----------|
| Identical sequences | distance = 0 (all methods) |
| All-gap alignment / empty sequences | distance = **0** — no comparable sites *and* no differences (`0/n → 0`); read as "no evidence of divergence", not an error |
| Unequal-length sequences | throw `ArgumentException` (aligned sequences are a precondition) |
| Null sequences | throw `ArgumentNullException` (API contract) |

## Worked oracles (from the source datasets)

- **Single difference** — `ACGTACGT` vs `TCGTACGT` (1 diff at pos 0): **Hamming = 1**, **p = 0.125**,
  **JC69 ≈ 0.137** (`-0.75·ln(1 − 4·0.125/3)`).
- **JC69 known values** — `p=0 → d=0`; `p=0.25 → d = -0.75·ln(2/3) ≈ 0.304`.
- **Pure transition** — `ACGT` vs `GCGT` (A→G at pos 0): `S=1/4, V=0` → **K2P = -0.5·ln(0.5) ≈ 0.34657**.
- **Pure transversion** — `ACGT` vs `CCGT` (A→C at pos 0): `S=0, V=1/4` →
  **K2P = -0.5·ln(0.75·√0.5) ≈ 0.31713** (note transversion < transition distance at equal p, the K2P
  signature).
- **Mixed** — `ACGTACGT` vs `GCGTTCGT` (1 transition + 1 transversion): `S=1/8, V=1/8` →
  **K2P ≈ 0.30679**.
- **With a gap** — `ACG-ACGT` vs `ACGTACGT`: the gap column is dropped, leaving **7 comparable sites**.

## Relationship to the rest of the PHYLO family

This unit is the **upstream substrate** of the family, not a tree method itself: it emits the distance
matrix that the [[distance-based-tree-construction]] step (PHYLO-TREE-001, **UPGMA/NJ**) turns into a
`PhyloNode` tree. [[phylogenetic-bootstrap-support]] resamples the
*alignment columns* and re-runs this distance step + tree builder per replicate (bootstrap **wraps**
the distance calculation); [[tree-comparison-metrics]] operates on the finished trees. It is **distinct
from** [[phylogenetic-marker-selection]] (PANGEN-MARKER-001), which selects the informative *columns* a
distance matrix would be computed over but computes no distance, and from
[[tumor-phylogeny-clonal-tree-reconstruction]] (ONCO-PHYLO-001), the oncology CCF-constraint tree
builder that computes **no** distance matrix at all. It is separate again from the k-mer / edit /
alignment-based sequence-similarity distances used elsewhere in the library (those are alignment-free or
scoring-based; these four are **substitution-model** evolutionary distances over an existing alignment).

## Documented assumptions (source-consistent)

1. **Empty / all-gap input → distance 0** (rather than NaN or throw) — the mathematical limit `0/n → 0`
   with zero differences; a defined API-contract behavior for "no comparable sites, no divergence".
2. **Pairwise deletion for gaps and ambiguity codes** — gaps and non-ACGT IUPAC bases are skipped
   per-position, so `comparableSites` can be fewer than the alignment length. (Complete-deletion and
   gap-as-mismatch are literature alternatives this unit does not take.)

No source contradictions: the JC69 / K2P formulas, the symmetry / zero-diagonal / correction-ordering
properties, and the saturation limits are the standard textbook definitions.

## Reference tools

Definitions trace to **Jukes & Cantor (1969)** *"Evolution of Protein Molecules"* (in Munro (ed.),
*Mammalian Protein Metabolism*, pp. 21–132 — the JC69 correction), **Kimura (1980)** *"A simple method
for estimating evolutionary rates of base substitutions…"* (J. Mol. Evol. 16:111–120 — the K2P
transition/transversion model), **Felsenstein (2004)** *Inferring Phylogenies* (Sinauer), and the
Wikipedia articles **Models of DNA evolution**, **Substitution model**, and **Distance matrices in
phylogeny** (symmetric zero-diagonal matrix, time-reversibility, gaps ignored or counted as mismatches
by implementation choice). No source contradictions.
