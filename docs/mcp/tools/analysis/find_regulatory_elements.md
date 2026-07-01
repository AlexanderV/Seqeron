# find_regulatory_elements

Scan for built-in regulatory motifs in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_regulatory_elements` |
| **Method ID** | `MotifFinder.FindRegulatoryElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence for a built-in catalog of **regulatory motifs**: TATA box, CAAT
box, GC box (Sp1), prokaryotic −10 (Pribnow) and −35 boxes, Kozak, Shine-Dalgarno,
poly(A) signal, E-box, AP-1, NF-κB and CREB sites. Matching uses IUPAC-degenerate
patterns, so ambiguity codes in a motif match the corresponding base sets. Each hit
reports the motif name, 0-based position, matched subsequence, pattern and a
description.

## Core Documentation Reference

- Source: [MotifFinder.cs#L685](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L685)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ name, position, sequence, pattern, description }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: TATA box

**User Prompt:**
> Find regulatory elements in "GGTATAAAGG".

**Expected Tool Call:**
```json
{
  "tool": "find_regulatory_elements",
  "arguments": { "sequence": "GGTATAAAGG" }
}
```

**Response:**
```json
{ "items": [ { "name": "TATA Box", "position": 2, "sequence": "TATAAA", "pattern": "TATAAA", "description": "Eukaryotic core promoter element" } ] }
```
The TATA box consensus TATAAA is found at position 2.

### Example 2: No regulatory element

**User Prompt:**
> Regulatory elements in "AAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "find_regulatory_elements",
  "arguments": { "sequence": "AAAAAAAA" }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · P) for P catalog patterns.
- **Space Complexity:** O(number of hits).

## See Also

- [find_known_motifs](find_known_motifs.md) — search a user-provided motif set
- [find_degenerate_motif](find_degenerate_motif.md) — single IUPAC pattern search
