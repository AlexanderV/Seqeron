# predict_rna_structure

Predict RNA secondary structure by greedy stem-loop selection.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_rna_structure` |
| **Method ID** | `RnaSecondaryStructure.PredictStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts an RNA **secondary structure** by greedily selecting non-overlapping
stem-loops (most stable first). Returns the echoed sequence, dot-bracket notation, the
base pairs, the selected stem-loops, any pseudoknots, and the total minimum free
energy. Structural elements do not overlap (nested pairs only).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L1839](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L1839)

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
| `minimumFreeEnergy` | number | Total MFE (kcal/mol) |

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
{ "sequence": "GGGAAAACCC", "dotBracket": "(((....)))" }
```
A single 3bp hairpin around the AAAA loop.

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

- **Time Complexity:** O(n²) stem-loop enumeration + greedy selection.
- **Space Complexity:** O(number of stem-loops).

## See Also

- [find_stem_loops](find_stem_loops.md) — the candidate hairpins
- [minimum_free_energy](minimum_free_energy.md) — Zuker MFE
- [parse_dot_bracket](parse_dot_bracket.md) — extract pairs from the notation
