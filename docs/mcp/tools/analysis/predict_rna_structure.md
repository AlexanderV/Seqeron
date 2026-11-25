# predict_rna_structure

Assemble an RNA secondary structure from hairpin stem-loops.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_rna_structure` |
| **Method ID** | `RnaSecondaryStructure.PredictStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assembles an RNA **secondary structure** from the hairpin stem-loops of `find_stem_loops`
(perfect stems of ≥ `minStemLength` pairs, loops of `minLoopSize`..`maxLoopSize` nt): the
lowest-free-energy set of non-overlapping stem-loops, chosen exactly by the Zuker–Stiegler
exterior-loop recursion over the candidates (each scored in context, d2 exterior terms included).
Stem-loops that would be less stable than the open chain are never selected. Returns the
normalised sequence, dot-bracket, base pairs, the selected stem-loops, pseudoknots (always
empty) and `minimumFreeEnergy` = the free energy of the returned structure (= `rna_structure_energy`,
≤ 0). The search space has no multiloops, bulges or interior loops, so the result can be less stable
than the true MFE fold — use [rna_mfe_structure](rna_mfe_structure.md) for that.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) (`PredictStructure`)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (min length 1) |
| `minStemLength` | integer | No | Minimum stem length (default 3) |
| `minLoopSize` | integer | No | Minimum loop size (default 3) |
| `maxLoopSize` | integer | No | Maximum loop size (default 10) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Echoed input sequence |
| `dotBracket` | string | Dot-bracket structure |
| `basePairs` | array | Base pairs |
| `stemLoops` | array | Selected stem-loops |
| `pseudoknots` | array | Detected pseudoknots |
| `minimumFreeEnergy` | number | Free energy of the returned structure (kcal/mol, ≤ 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Simple hairpin

**User Prompt:**
> Predict the structure of "GGGAAAACCC".

**Expected Tool Call:**
```json
{
  "tool": "predict_rna_structure",
  "arguments": { "rnaSequence": "GGGAAAACCC", "minStemLength": 3, "minLoopSize": 4, "maxLoopSize": 4 }
}
```

**Response:**
```json
{ "sequence": "GGGAAAACCC", "dotBracket": "(((....)))", "minimumFreeEnergy": -2.02 }
```
A single 3bp hairpin around the AAAA loop (ViennaRNA `eval_structure` −2.02).

### Example 2: Unstructured

**User Prompt:**
> Predict the structure of "AAAAAAAAAA".

**Expected Tool Call:**
```json
{
  "tool": "predict_rna_structure",
  "arguments": { "rnaSequence": "AAAAAAAAAA" }
}
```

**Response:**
```json
{ "dotBracket": ".........." }
```
No stable pairs ⇒ all-dot structure.

## Performance

- **Time Complexity:** O(n·(maxLoopSize − minLoopSize + 1)) candidates, each scored in O(n) → O(n²); selection O(n + candidates).
- **Space Complexity:** O(n + number of candidates).

## See Also

- [find_stem_loops](find_stem_loops.md) — the candidate hairpins
- [minimum_free_energy](minimum_free_energy.md) — Zuker MFE
- [rna_mfe_structure](rna_mfe_structure.md) — MFE structure
- [parse_dot_bracket](parse_dot_bracket.md) — extract pairs from the notation
