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
- **AvoidRareCodons** — only replace codons below the rare-codon threshold (with the most frequent synonym). The historical misspelling **AvoidRareCodeons** is the same value and stays accepted; strategy names are matched case-insensitively (review 2026-09, F55).
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
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard, identical to omitting it; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the encoded protein and the synonymous-codon families (DNA Chisel `EnforceTranslation(genetic_table=…)`); under table 2 AGA/AGG are stops and never used for Arg, UGA is Trp, AUA is Met. One of 1–6, 9–16, 21–33. |
| `avoid_sites` | array \| null | No | `null`/omitted (default): no site avoidance. A list of restriction recognition sequences (DNA or RNA, IUPAC codes allowed — same form as [remove_restriction_sites](remove_restriction_sites.md), which also takes sequences rather than enzyme names) runs the site-removal pass **after** codon selection, the GC pass and the structure pass, so no later pass can re-introduce a site (DNA Chisel resolves the hard `AvoidPattern` constraints last). Both strands are cleared; each removal changes one codon; an occurrence no single synonymous substitution can remove is left in place. For the strategies with a GC target (`BalancedOptimization`, `MinimizeSecondary`) the pass is **GC-constraint-aware** (review 2026-09, F60): it keeps `gc_target_min`/`gc_target_max`, and every `gc_window` window, satisfied whenever a site-removing synonymous codon can — ranking candidates by GC feasibility (or, when the input already breaches, by not increasing the breach), then codon frequency, then the leftmost codon, then NCBI codon order. If no candidate can, the least-breach candidate is applied and the site is still removed (site avoidance stays the hard constraint; dnachisel would raise `NoSolutionError`). The reported CAI, GC fractions and codon changes always describe the final sequence (review 2026-09, F54/F60). |

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
| 1006 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| 1007 | Restriction site '…' contains a character that is not an IUPAC nucleotide code (from `avoid_sites`) |

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

### Example 3: Vertebrate mitochondrial gene (`genetic_code` 2)

`AUG CGU AGA CGC AGG UGA AUA UAA` with `strategy` MaximizeCAI and the E. coli K12 preset: under the default Standard code the result is `AUGCGCCGCCGCCGCUGAAUUUAA` (AGA/AGG recoded as Arg CGC, AUA → AUU Ile); with `genetic_code: 2` it is `AUGCGCAGACGCAGGUGGAUGUAA` — AGA/AGG are table-2 stops and kept, UGA is Trp (→ UGG), AUA is Met (→ AUG); `proteinSequence` = `MR*R*WM*`, the table-2 protein of the input (DNA Chisel `EnforceTranslation(genetic_table="Vertebrate Mitochondrial")` families; locked in `CodonOptimizer_GeneticCode_Tests`).

### Example 4: Chained restriction-site avoidance (`avoid_sites`)

`optimize_codons("ATGCTACAACACATGACAATGGCAAAGTAA", { preset: "EColiK12" }, strategy: "MaximizeCAI", avoid_sites: ["CTGCAG", "CATATG", "CCATGG"])` — MaximizeCAI alone would *create* all three sites (Leu CUG + Gln CAG = PstI, His CAU + Met AUG = NdeI, Thr ACC + Met AUG = NcoI). With `avoid_sites` the result is `AUGCUGCAACACAUGACGAUGGCGAAAUAA` (protein `MLQHMTMAK*`, `optimizedCAI` 0.793311 instead of 1.0, 4 changed codons) — byte-identical to the sequence dnachisel 3.2.16 returns for `DnaOptimizationProblem(constraints=[AvoidPattern(...)×3, EnforceTranslation()], objectives=[CodonOptimize(species='e_coli', method='use_best_codon')])`, whose `AvoidPattern(...).evaluate` scores our output 0 (pass) for every site.

### Example 5: Empty input

`optimize_codons("", { preset: "EColiK12" })` throws `ArgumentException`.

## See Also

- [cai_from_organism_table](cai_from_organism_table.md), [find_rare_codons](find_rare_codons.md), [build_codon_table](build_codon_table.md)
