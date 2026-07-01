# fit_heaps_law

Fit Heaps' law to the pan-genome new-gene-discovery curve.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `fit_heaps_law` |
| **Method ID** | `PanGenomeAnalyzer.FitHeapsLaw` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Clusters the genomes into ortholog groups, builds the presence/absence matrix, then fits
Heaps' law `n(N) = Intercept · N^(-Alpha)` to the number of *new* gene clusters contributed by
the N-th genome, pooled over `permutations` random genome orderings (micropan `heaps()`,
Tettelin et al. 2008). The pan-genome is **open** when `Alpha < 1`. Permutations use a fixed
seed, and the first permutation is the natural input order, so `permutations = 1` is exactly
reproducible.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L588](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L588)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `identityThreshold` | number | No | Ortholog-clustering identity (default 0.9). |
| `permutations` | integer | No | Random orderings to pool over (default 100). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `intercept` | number | Heaps' law intercept K ∈ [0, 10000]. |
| `alpha` | number | Decay exponent α ∈ [0, 2]. |
| `isOpen` | boolean | `true` when α < 1 (open pan-genome). |

## Errors

None. Fewer than 2 genomes yields a degenerate closed fit `(0, 0, false)`.

## Example

**User Prompt:**
> Is this pan-genome open? Each genome adds exactly one new gene.

With `permutations = 1` on a constant new-gene curve:

**Response:**
```json
{ "intercept": 1.0, "alpha": 0.0, "isOpen": true }
```

A flat new-gene curve fits α = 0, K = 1 → open.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `FitHeapsLaw`
- Tettelin H. et al. (2008) Curr. Opin. Microbiol. 11:472; micropan `heaps()`.

## See Also

- [construct_pangenome](construct_pangenome.md)
- [gene_presence_absence_matrix](gene_presence_absence_matrix.md)
