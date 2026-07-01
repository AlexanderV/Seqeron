# align_mirna_to_target

Align a miRNA against a target sequence and compute duplex statistics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `align_mirna_to_target` |
| **Method ID** | `MiRnaAnalyzer.AlignMiRnaToTarget` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Aligns a miRNA to a target sequence in the biologically correct antiparallel orientation:
miRNA position `i` (read 5'→3') is paired against the target base read 3'→5', i.e.
`target[len-1-i]`. Each position is scored as a Watson-Crick match (`|`), a G-U wobble
(`:`, Crick 1966), or a mismatch (` `). The duplex free energy is estimated by summing
Turner 2004 nearest-neighbor stacking energies over consecutive paired positions. Inputs are
case-insensitive and DNA `T` is treated as RNA `U`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L364](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L364)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence (min length: 1) |
| `targetSequence` | string | Yes | Target nucleotide sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `miRnaSequence` | string | Normalised (uppercase RNA) miRNA sequence |
| `targetSequence` | string | Normalised (uppercase RNA) target sequence |
| `alignmentString` | string | Per-position symbols: `|` Watson-Crick, `:` G-U wobble, ` ` mismatch |
| `matches` | integer | Number of Watson-Crick paired positions |
| `mismatches` | integer | Number of unpaired positions |
| `guWobbles` | integer | Number of G-U wobble paired positions |
| `gaps` | integer | Number of gaps (always 0; ungapped alignment) |
| `freeEnergy` | number | Turner 2004 nearest-neighbor duplex free energy (kcal/mol, 37 °C) |

## Errors

| Code | Message |
|------|---------|
| 1001 | miRNA sequence cannot be null or empty |
| 1001 | Target sequence cannot be null or empty |

## Examples

### Example 1: Perfect antiparallel complement

**User Prompt:**
> Align miRNA "AAAA" to target "UUUU".

**Expected Tool Call:**
```json
{
  "tool": "align_mirna_to_target",
  "arguments": { "miRnaSequence": "AAAA", "targetSequence": "UUUU" }
}
```

**Response:**
```json
{
  "matches": 4,
  "mismatches": 0,
  "guWobbles": 0,
  "gaps": 0,
  "alignmentString": "||||"
}
```

### Example 2: G-U wobble pairing

**User Prompt:**
> Align miRNA "GGGG" to target "UUUU".

**Response:**
```json
{
  "matches": 0,
  "mismatches": 0,
  "guWobbles": 4,
  "gaps": 0,
  "alignmentString": "::::"
}
```

## Performance

- **Time Complexity:** O(n) where n is min(miRNA length, target length)
- **Space Complexity:** O(n)

## See Also

- [find_mirna_target_sites](find_mirna_target_sites.md) - Scan an mRNA for miRNA target sites
- [can_pair](can_pair.md) - Whether two bases can pair
- [is_wobble_pair](is_wobble_pair.md) - Whether a pair is a G-U wobble
