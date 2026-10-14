# dinucleotide_ratios

Observed/expected dinucleotide relative abundance (odds ratios).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dinucleotide_ratios` |
| **Method ID** | `SequenceStatistics.CalculateDinucleotideRatios` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the dinucleotide relative abundance (Karlin odds ratio) `ρ_XY = f_XY / (f_X · f_Y)`, where
`f_XY` is the dinucleotide frequency and `f_X`, `f_Y` are single-base frequencies. `ρ = 1`
means no bias; `ρ > 1` over-representation, `ρ < 1` under-representation (Karlin & Burge 1995).

- **Single strand** (default, `strandSymmetric = false`): alphabet A/C/G/T/U; other symbols (N,
  IUPAC codes, gaps) are excluded from counts and denominators. Identical to seqinr `rho(seq)` and,
  for pure ACGT input, to EMBOSS `compseq -word 2 -calcfreq` Obs/Exp.
- **Strand-symmetrized** (`strandSymmetric = true`): Karlin's genomic-signature `ρ*_XY`, computed from the
  sequence together with its reverse complement (`f*_A = f*_T = (f_A + f_T)/2`, `f*_GT = (f_GT + f_AC)/2`, …;
  no junction dinucleotide). Alphabet A/C/G/T (U excluded like N).
- **Representation classes** (`includeRepresentation = true`): adds Karlin's class per returned dinucleotide —
  `UnderRepresented` (ρ ≤ 0.78), `Normal`, `OverRepresented` (ρ ≥ 1.23) (Karlin & Burge 1995; Karlin 1998).

Only dinucleotides that occur (on either strand for ρ*) are keys; an absent key means ρ = 0 when both
bases occur and undefined otherwise. Sequences with fewer than two counted bases yield an empty map.

This odds ratio normalises by the counted dinucleotides (≈ N − 1). It is **not** the Gardiner-Garden &
Frommer CpG observed/expected ratio used for CpG-island calling (normalised by the window length N,
threshold 0.6) — that lives in `EpigeneticsAnalyzer` (annotation server CpG-island tools).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L1308](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L1308)
- Classes: `SequenceStatistics.ClassifyDinucleotideRatios` (Karlin thresholds 0.78 / 1.23)
- Evidence: `docs/Evidence/SEQ-DINUC-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1, case-insensitive) |
| `strandSymmetric` | boolean | No | `true`: Karlin ρ* over the sequence and its reverse complement; `false` (default): single-strand ρ |
| `includeRepresentation` | boolean | No | `true`: also return the Karlin class per dinucleotide (default `false`, field omitted) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ratios` | object | Map of dinucleotide → odds ratio ρ (or ρ*) |
| `representation` | object | Only with `includeRepresentation`: dinucleotide → `UnderRepresented` / `Normal` / `OverRepresented` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Alternating AT

**Input:** `{ "sequence": "ATAT" }`
→ f_A = f_T = 0.5; f_AT = 2/3, f_TA = 1/3 →
`ρ_AT = (2/3)/0.25 = 2.6667`, `ρ_TA = (1/3)/0.25 = 1.3333`.

**Response:**
```json
{ "ratios": { "AT": 2.6667, "TA": 1.3333 } }
```

### Example 2: Homopolymer

**Input:** `{ "sequence": "AAAA" }`
→ f_A = 1, f_AA = 1 → `ρ_AA = 1 / (1·1) = 1.0`.

### Example 3: Karlin strand-symmetrized ρ* with classes

**Input:** `{ "sequence": "ATGCGCGT", "strandSymmetric": true, "includeRepresentation": true }`
→ seqinr `rho` of the sequence + separator + reverse complement (exact fractions): AT = 256/63,
CG = GC = 512/175, AC = CA = GT = TG = 128/105.

**Response:**
```json
{
  "ratios": { "AC": 1.2190, "AT": 4.0635, "CA": 1.2190, "CG": 2.9257, "GC": 2.9257, "GT": 1.2190, "TG": 1.2190 },
  "representation": { "AC": "Normal", "AT": "OverRepresented", "CA": "Normal", "CG": "OverRepresented", "GC": "OverRepresented", "GT": "Normal", "TG": "Normal" }
}
```

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct dinucleotides).

## See Also

- [dinucleotide_frequencies](dinucleotide_frequencies.md)
