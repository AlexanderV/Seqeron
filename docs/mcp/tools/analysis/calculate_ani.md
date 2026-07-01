# calculate_ani

Average Nucleotide Identity (ANI) between two genomes.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `calculate_ani` |
| **Method ID** | `ComparativeGenomics.CalculateANI` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the Average Nucleotide Identity (ANI) between two genome sequences using
the ANIb definition of Goris et al. (2007). The first genome is cut into consecutive
non-overlapping fragments of `fragmentSize` nt; each fragment is placed by its best
ungapped local match against genome 2. A fragment contributes only if its identity
(recalculated over the whole fragment length) exceeds `minFragmentIdentity` **and**
its alignable region covers at least 70% of the fragment. ANI is the mean identity
of the qualifying fragments, returned as a fraction in `[0, 1]`. ANI ≈ 0.95
corresponds to the 70% DDH species boundary. Returns 0 when either sequence is empty,
no full fragment fits, or no fragment qualifies.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L1076](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L1076)
- Evidence: `docs/Evidence/COMPGEN-ANI-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Sequence` | string | Yes | Query genome (fragmented) — min length 1 |
| `genome2Sequence` | string | Yes | Reference genome — min length 1 |
| `fragmentSize` | integer | No | Fragment length in nt (default 1000, > 0) |
| `minFragmentIdentity` | number | No | Minimum per-fragment identity to keep a match (default 0.7) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ani` | number | Average Nucleotide Identity as a fraction in [0, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | Genome sequence cannot be null or empty |
| 1003 | Fragment size must be positive |

## Examples

### Example 1: Identical genomes

**User Prompt:**
> ANI of "AAAACCCCGGGGTTTT" against itself with 4-nt fragments.

**Expected Tool Call:**
```json
{
  "tool": "calculate_ani",
  "arguments": {
    "genome1Sequence": "AAAACCCCGGGGTTTT",
    "genome2Sequence": "AAAACCCCGGGGTTTT",
    "fragmentSize": 4
  }
}
```

**Response:**
```json
{ "ani": 1.0 }
```
Every 4-nt fragment is a perfect substring ⇒ mean identity 1.0.

### Example 2: One substitution

**User Prompt:**
> ANI of "AAAACCCCGGGGTTTA" vs "AAAACCCCGGGGTTTT", 4-nt fragments.

**Expected Tool Call:**
```json
{
  "tool": "calculate_ani",
  "arguments": {
    "genome1Sequence": "AAAACCCCGGGGTTTA",
    "genome2Sequence": "AAAACCCCGGGGTTTT",
    "fragmentSize": 4
  }
}
```

**Response:**
```json
{ "ani": 0.9375 }
```
Last fragment `TTTA` matches `TTTT` 3/4 = 0.75; ANI = (1+1+1+0.75)/4 = 0.9375.

## Performance

- **Time Complexity:** O(|query|/fragmentSize · |reference|) for ungapped placement.
- **Space Complexity:** O(qualifying fragments).

## See Also

- [compare_genomes](compare_genomes.md)
- [find_orthologs](find_orthologs.md)
