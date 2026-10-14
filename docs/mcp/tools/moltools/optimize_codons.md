# optimize_codons

Optimize a coding sequence's codons for expression in a target organism.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `optimize_codons` |
| **Method ID** | `CodonOptimizer.OptimizeSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Rewrites the synonymous codons of a coding sequence for a target organism (an `EColiK12`/`Yeast`/`Human` preset or an inline custom RNA-alphabet table) under one of five strategies. The input is upper-cased, converted T→U, and trimmed to whole codons; stop codons and single-codon amino acids (Met, Trp) are never changed. Returns the original and optimized RNA, the translated protein, original/optimized CAI (Sharp & Li 1987), original/optimized GC fraction, the number of changed codons, and each `(position, original, optimized)` change.

- **MaximizeCAI** — pick the most frequent synonymous codon; ties break by NCBI codon order (identical to DNA Chisel `use_best_codon`).
- **BalancedOptimization** (default) — as MaximizeCAI, then GC-balance into `[gc_target_min, gc_target_max]` with further synonymous swaps. The GC pass constrains the **global** GC fraction by default (DNA Chisel `EnforceGCContent(mini, maxi)`); with `gc_window = w` it instead constrains **every sliding window** of `w` nucleotides (DNA Chisel `EnforceGCContent(mini, maxi, window=w)`), resolved by a deterministic greedy synonymous-codon repair (review 2026-09, F43).
- **AvoidRareCodeons** — only replace codons below the rare-codon threshold (with the most frequent synonym).
- **MinimizeSecondary** — balanced selection and GC pass, **followed by** the sliding-window structure pass of `reduce_secondary_structure`, scored by `secondary_scoring`: `PairCount` (default, the canonical-pair fraction heuristic, ~0.06 s per kilobase) or `Turner2004Mfe` (Turner 2004 nearest-neighbour free energy, the ViennaRNA RNAfold energy model, ~9 s per kilobase).
- **HarmonizeExpression** — match the target codon-usage profile as closely as integer codon counts allow (DNA Chisel `match_codon_usage`); **deterministic** since review 2026-09.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L469](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L469)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12`/`Yeast`/`Human`) or inline custom table. |
| `strategy` | enum | No | Optimization strategy (default `BalancedOptimization`). |
| `gc_target_min` | number | No | Lower GC fraction bound in [0,1] (default 0.40). |
| `gc_target_max` | number | No | Upper GC fraction bound in [0,1] (default 0.60). |
| `rare_codon_threshold` | number | No | Rare-codon frequency threshold (default 0.15). |
| `secondary_scoring` | enum | No | Scoring of the MinimizeSecondary structure pass: `PairCount` (default) or `Turner2004Mfe`. |
| `gc_window` | integer \| null | No | `null` (default): global GC constraint only. A positive window size in nucleotides constrains every sliding window of that size to `[gc_target_min, gc_target_max]`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `originalSequence` / `optimizedSequence` | string | RNA sequences (T→U). |
| `proteinSequence` | string | Translated protein (`*` = stop). |
| `originalCAI` / `optimizedCAI` | number | Codon Adaptation Index (0–1). |
| `gcContentOriginal` / `gcContentOptimized` | number | GC fraction. |
| `changedCodons` | integer | Number of codons changed. |
| `changes` | array | `(position, original, optimized)` codon changes. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | GC target minimum/maximum must be a fraction in [0, 1] |
| 1003 | GC target minimum must not exceed the maximum |
| 1004 | Unknown codon-usage preset |
| 1005 | GC window must be a positive number of nucleotides |

## Examples

### Example 1: MaximizeCAI on rare E. coli leucines

`ATG TTA CTT CTA TAA` → all three rare Leu codons become `CUG` (E. coli's preferred, freq 0.50), optimized CAI = 1.0.

**Input:** `{ "coding_sequence": "ATGTTACTTCTATAA", "target_organism": { "preset": "EColiK12" }, "strategy": "MaximizeCAI" }`

**Response (abridged):**
```json
{
  "originalSequence": "AUGUUACUUCUAUAA",
  "optimizedSequence": "AUGCUGCUGCUGUAA",
  "proteinSequence": "MLLL*",
  "optimizedCAI": 1.0,
  "changedCodons": 3,
  "changes": [
    { "position": 3, "original": "UUA", "optimized": "CUG" },
    { "position": 6, "original": "CUU", "optimized": "CUG" },
    { "position": 9, "original": "CUA", "optimized": "CUG" }
  ]
}
```

### Example 2: Windowed GC constraint

`optimize_codons("ATGAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTTAA", { preset: "EColiK12" }, gc_window: 20)` keeps every 20-nt window inside `[0.40, 0.60]` — DNA Chisel `EnforceGCContent(0.40, 0.60, window=20).evaluate(...)` scores the result 0 (the globally balanced default output still scores below 0 on 20-nt windows).

### Example 3: Empty input

`optimize_codons("", { preset: "EColiK12" })` throws `ArgumentException`.

## See Also

- [cai_from_organism_table](cai_from_organism_table.md), [find_rare_codons](find_rare_codons.md), [build_codon_table](build_codon_table.md)
