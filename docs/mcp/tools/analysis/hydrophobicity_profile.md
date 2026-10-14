# hydrophobicity_profile

Sliding-window Kyte-Doolittle hydropathy profile of a protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `hydrophobicity_profile` |
| **Method ID** | `SequenceStatistics.CalculateHydrophobicityProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **Kyte-Doolittle hydropathy profile** (Kyte & Doolittle, 1982) as ExPASy ProtScale /
Biopython `ProteinAnalysis.protein_scale(kd, window, edge)`: the weighted mean hydropathy over each
sliding window of `windowSize` residues. Positive values indicate hydrophobic stretches (candidate
membrane-spanning or buried regions), negative values hydrophilic ones. The profile has
`N − windowSize + 1` values; empty when the window exceeds the sequence length. Input is case-insensitive.

- **`edgeWeight`** (default 1 = unweighted mean): ProtScale "relative weight of the window edges" — the
  window weights rise linearly from `edgeWeight` at both edges to 1 at the centre (Biopython `_weight_list`).
  Must be in [0, 1]; < 1 requires an odd window in the Standard convention.
- **`convention`**: `Standard` (default) — each position has its own weight, a residue without a value
  contributes 0 and keeps its weight; `Biopython` — line-for-line Biopython 1.88 `protein_scale`: for an even
  window residue W/2 is counted twice (divisor W + 1 at edge 1), a symmetric pair containing an unknown
  residue is dropped as a whole; window ≥ 2.
- **`strict`**: reject any symbol without a scale value (Biopython `gravy()` KeyError behaviour).
- **`normalizeScale`**: min–max normalise the Kyte-Doolittle scale to 0..1 before profiling (ProtScale
  "normalize the scale from 0 to 1").

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L676](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L676) (defaults / `edgeWeight`), [#L753](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L753) (`convention` / `strict` / `normalizeScale`)
- Algorithm: `docs/algorithms/Statistics/Hydrophobicity_Analysis.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 9, ≥ 1; ≥ 2 with `Biopython`) |
| `edgeWeight` | number | No | Relative weight of the window edges in [0, 1] (default 1) |
| `convention` | string | No | `Standard` (default) or `Biopython` (case-insensitive) |
| `strict` | boolean | No | Reject symbols without a scale value (default `false`) |
| `normalizeScale` | boolean | No | Normalise the scale to 0..1 first (default `false`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | array of number | Weighted mean hydropathy per window (length `N − windowSize + 1`) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |
| 1002 | Edge weight must be in [0, 1]. |
| 1003 | An edge-weighted window must have an odd size (a central residue). |
| 1003 | Unknown window convention '…'. Expected one of: Standard, Biopython. |
| 1003 | Symbol '…' has no value in the scale (strict mode). |

## Examples

### Example 1: Hydrophobic poly-isoleucine (window 3)

**User Prompt:**
> Hydrophobicity profile of "IIIII" with window 3.

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity_profile",
  "arguments": { "proteinSequence": "IIIII", "windowSize": 3 }
}
```

**Response:**
```json
{ "values": [4.5, 4.5, 4.5] }
```
Isoleucine has KD value 4.5, so every window averages 4.5.

### Example 2: Mixed residues (window 2)

**User Prompt:**
> Hydrophobicity of "AV" with window 2.

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity_profile",
  "arguments": { "proteinSequence": "AV", "windowSize": 2 }
}
```

**Response:**
```json
{ "values": [3.0] }
```
(A = 1.8, V = 4.2) ⇒ mean = 3.0.

### Example 3: Edge-weighted window (ProtScale / Biopython)

**Input:** `{ "proteinSequence": "FLIVAG", "windowSize": 5, "edgeWeight": 0.4 }`

**Response:** `{ "values": [3.73125, 3.1156250000000005] }`
(Biopython 1.88 `ProteinAnalysis("FLIVAG").protein_scale(kd, 5, 0.4)` = 3.7312499999999997, 3.1156249999999996 —
equal up to floating-point summation order; weights 0.4, 0.7, 1, 0.7, 0.4.)

### Example 4: Biopython window convention with an unknown residue

**Input:** `{ "proteinSequence": "AXCDE", "windowSize": 5, "convention": "Biopython" }` → `{ "values": [0.16] }`
(Biopython drops the pair X/D; the Standard convention gives −0.54.)

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(n) for the profile.

## See Also

- [predict_transmembrane_helices](predict_transmembrane_helices.md) — hydropathy-based TM prediction
- [predict_chou_fasman](predict_chou_fasman.md) — secondary-structure propensities
