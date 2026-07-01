# find_pre_mirna_hairpins

Identify pre-miRNA hairpin candidates in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_pre_mirna_hairpins` |
| **Method ID** | `MiRnaAnalyzer.FindPreMiRnaHairpins` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence for stem-loop (hairpin) windows consistent with pre-miRNA structure — sufficient stem
pairing, a loop of 3–25 nt, and total length within `[minHairpinLength, maxHairpinLength]`. Each candidate
reports its span, the extracted mature/star arm sequences, a dot-bracket structure, and a Turner (2004)
nearest-neighbour free-energy estimate. This is a simplified consecutive-pairing model; real pre-miRNAs
with internal bulges may not be detected.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L1744](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L1744)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence to scan |
| `minHairpinLength` | integer | No | Minimum hairpin length in nt (default 55) |
| `maxHairpinLength` | integer | No | Maximum hairpin length in nt (default 120) |
| `matureLength` | integer | No | Mature miRNA length in nt (default 22) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hairpins[].start` | integer | Hairpin start |
| `hairpins[].end` | integer | Hairpin end |
| `hairpins[].sequence` | string | Full hairpin sequence |
| `hairpins[].matureSequence` | string | Extracted mature arm |
| `hairpins[].starSequence` | string | Extracted star arm |
| `hairpins[].structure` | string | Dot-bracket secondary structure |
| `hairpins[].freeEnergy` | number | Turner 2004 nearest-neighbour free energy |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: A 57 nt hairpin

A 25 bp-stem / 7 nt-loop / 25 bp-stem sequence with `minHairpinLength = 55` yields two candidate windows
(the full 57 nt window at 0–56 and a 55 nt sub-window at offset 1):

**Response:**
```json
{ "hairpins": [ { "start": 0, "end": 56 }, { "start": 1, "end": 56 } ] }
```

### Example 2: Below minimum length

A 47 nt hairpin is rejected against a 55 nt minimum:

**Response:**
```json
{ "hairpins": [] }
```

## Performance

- **Time Complexity:** O(n · L) for sequence length n and window length L
- **Space Complexity:** O(L)

## See Also

- [create_mirna](create_mirna.md) — build a mature miRNA record
- [align_mirna_to_target](align_mirna_to_target.md) — miRNA-target duplex
