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

Scans a 200 bp window upstream of `target_start` for forward-primer candidates and a 200 bp window downstream of `target_end` for reverse-primer candidates (evaluated on the reverse complement), keeps only those that pass every quality gate in `PrimerParameters`, and returns the highest-scoring valid candidate on each side. The pair is flagged compatible when the two Tm values differ by ≤ 5 °C and the primers do not form a 3′ primer-dimer. `product_size = reverse.Position + reverse.Length − forward.Position`.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L40](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L40)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | DNA template (A/C/G/T), non-empty. |
| `target_start` | integer | Yes | 0-based inclusive start of the target region (≥ 0). |
| `target_end` | integer | Yes | 0-based inclusive end of the target region (`target_start < target_end < template.Length`). |
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

Forward region of `GAACTCGT` units, a 50 bp poly-T target, and a reverse region of `TCCGAAGT` units.

**Input:** `{ "target_start": 100, "target_end": 150 }`

**Response (abridged):**
```json
{
  "isValid": true,
  "productSize": 180,
  "forward": { "position": 0, "length": 25, "gcContent": 52.0, "meltingTemperature": 59.3 },
  "reverse": { "position": 155, "length": 25, "gcContent": 52.0, "meltingTemperature": 59.3 }
}
```

### Example 2: Invalid target region

`design_primers("ACGT…", 150, 100)` (start ≥ end) throws `ArgumentException`.

## See Also

- [generate_primer_candidates](generate_primer_candidates.md), [evaluate_primer](evaluate_primer.md), [primer_dimer](primer_dimer.md)
