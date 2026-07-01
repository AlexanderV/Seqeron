# validate_probe

Validate a probe's specificity against reference sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `validate_probe` |
| **Method ID** | `ProbeDesigner.ValidateProbe` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans each reference sequence with ungapped `max_mismatches`-tolerant (Hamming) matching, counts the total hits, computes the probe's self-complementarity and a secondary-structure flag, and derives a `[0,1]` specificity score: **0 hits → 0.0**, **1 hit → 1.0**, **N hits → 1/N**. `IsValid` is true when there are no issues, or when hits ≤ 1 and self-complementarity ≤ 0.4. Matching is case-insensitive.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L857](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L857)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `probe_sequence` | string | Yes | Probe sequence to validate (non-null; empty → invalid). |
| `reference_sequences` | string[] | Yes | Reference sequences to scan (non-null). |
| `max_mismatches` | integer | No | Maximum allowed mismatches (≥ 0, default 3). |
| `self_complementarity_threshold` | number | No | Self-complementarity warning threshold (default 0.3). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isValid` | boolean | Overall validity verdict. |
| `specificityScore` | number | `[0,1]` specificity (0 hits → 0, 1 → 1, N → 1/N). |
| `offTargetHits` | integer | Total approximate-match hits across references. |
| `selfComplementarity` | number | Self-complementarity fraction `[0,1]`. |
| `hasSecondaryStructure` | boolean | Secondary-structure potential flag. |
| `issues` | string[] | Reported issues. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Probe sequence cannot be null |
| 1002 | Reference sequences cannot be null |
| 1003 | Maximum mismatches cannot be negative |

## Examples

### Example 1: Unique probe → specificity 1.0

Probe `ATCGATCGATCGATCGATCG` in `NNNNNATCGATCGATCGATCGATCGNNNN` with `max_mismatches = 0`: 1 hit → specificity 1.0.

**Input:** `{ "probe_sequence": "ATCGATCGATCGATCGATCG", "reference_sequences": ["NNNNNATCGATCGATCGATCGATCGNNNN"], "max_mismatches": 0 }`

**Response (abridged):** `{ "offTargetHits": 1, "specificityScore": 1.0 }`

### Example 2: Repetitive probe → 1/N specificity

10-mer poly-A in a 34-mer poly-A: 25 exact-match positions → specificity `1/25 = 0.04`.

## See Also

- [design_probes](design_probes.md), [design_tiling_probes](design_tiling_probes.md), [crispr_specificity_score](crispr_specificity_score.md)
