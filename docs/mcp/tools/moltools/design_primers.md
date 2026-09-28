# design_primers

Design a forward/reverse PCR primer pair flanking a target region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_primers` |
| **Method ID** | `PrimerDesigner.DesignPrimers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a 200 bp window upstream of `target_start` for forward-primer candidates and a 200 bp window downstream of `target_end` for reverse-primer candidates (evaluated on the reverse complement), keeps only those that pass every quality gate in `PrimerParameters` (Tm = Primer3-default SantaLucia 1998 NN Tm), and returns — as Primer3 does — the pair with the lowest pair penalty (sum of the Primer3 per-primer penalties |Tm − OptimalTm| + |length − OptimalLength|) among all pairs whose Tm values differ by ≤ 5 °C and that do not form a 3′ primer-dimer; ties are broken as Primer3's `compare_primer_pair`. When no pair is compatible, the individually best primers are returned with `isValid = false`. `product_size = reverse.Position + reverse.Length − forward.Position`. Cross-checked against primer3-py 2.3.1 `design_primers`.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L40](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L40)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | DNA template (A/C/G/T), non-empty. |
| `target_start` | integer | Yes | 0-based inclusive start of the target region (≥ 0). |
| `target_end` | integer | Yes | 0-based **exclusive** end of the target region (`target_start < target_end < template.Length`); primers never overlap `[target_start, target_end)`. |
| `parameters` | object | No | Optional `PrimerParameters`; defaults (18–25 nt, 40–60% GC, 57–63 °C Tm) when null. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `forward` | object \| null | Best forward `PrimerCandidate` (null if none valid). |
| `reverse` | object \| null | Best reverse `PrimerCandidate` (null if none valid). |
| `isValid` | boolean | True when a compatible pair was found. |
| `message` | string | Human-readable status. |
| `productSize` | integer | Amplicon size in bp. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Template cannot be null or empty |
| 1002 | Target start must be non-negative |
| 1003 | Target end must be within the template |
| 1004 | Target start must be strictly less than target end |

## Examples

### Example 1: Standard 258 bp template, target [100, 150)

Forward region of `GAACTCGT` units, a 50 bp poly-T target, and a reverse region of `TCCGAAGT` units. primer3-py returns the same pair (PRIMER_LEFT_0 = [76,20], PRIMER_RIGHT_0 = [173,20], PRIMER_PAIR_0_PENALTY = 0.6983).

**Input:** `{ "target_start": 100, "target_end": 150 }`

**Response (abridged):**
```json
{
  "isValid": true,
  "productSize": 98,
  "forward": { "sequence": "TCGTGAACTCGTGAACTCGT", "position": 76, "length": 20, "gcContent": 50.0, "meltingTemperature": 59.3 },
  "reverse": { "sequence": "CGGAACTTCGGAACTTCGGA", "position": 154, "length": 20, "gcContent": 55.0, "meltingTemperature": 60.0 }
}
```

### Example 2: Invalid target region

`design_primers("ACGT…", 150, 100)` (start ≥ end) throws `ArgumentException`.

## See Also

- [generate_primer_candidates](generate_primer_candidates.md), [evaluate_primer](evaluate_primer.md), [primer_dimer](primer_dimer.md)
