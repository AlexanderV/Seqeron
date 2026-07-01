# rna_reverse_complement

Reverse-complement of an RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `rna_reverse_complement` |
| **Method ID** | `MiRnaAnalyzer.GetReverseComplement` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the reverse complement over the RNA alphabet: `A↔U`, `G↔C`, DNA `T` is complemented as `A`, and
any unrecognised base maps to `N`. Input is case-insensitive; output is upper-case. This is the operation
used to derive a miRNA seed's target pattern (e.g. let-7a seed `GAGGUAG` → `CUACCUC`).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L297](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L297)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA (or DNA) nucleotide sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `reverseComplement` | string | Reverse complement (unknown bases → N) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: RNA reverse complement

**User Prompt:**
> Reverse-complement AUGC.

**Response:**
```json
{ "reverseComplement": "GCAU" }
```

### Example 2: Unknown base

```json
{ "tool": "rna_reverse_complement", "arguments": { "rnaSequence": "AXGC" } }
```

**Response:**
```json
{ "reverseComplement": "GCNU" }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n)

## See Also

- [can_pair](can_pair.md) — Watson-Crick / wobble pairing test
- [find_mirna_target_sites](find_mirna_target_sites.md) — seed-RC-based target search
