# find_motif_by_pattern

Find regex pattern matches in a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_motif_by_pattern` |
| **Method ID** | `ProteinMotifFinder.FindMotifByPattern` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds all **overlapping** matches of a .NET regular expression in a protein sequence
(case-insensitive). Each match reports its 0-based start/end, the matched
subsequence, the supplied motif name and pattern id, a heuristic score and an E-value.
Zero-width captures are skipped, and pathological (catastrophic-backtracking) patterns
yield no matches rather than hanging.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L198](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L198)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `regexPattern` | string | Yes | .NET regex pattern (min length 1) |
| `motifName` | string | No | Motif display name (default `Custom`) |
| `patternId` | string | No | Pattern identifier (default empty) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, sequence, motifName, pattern, score, eValue }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Regex pattern cannot be null or empty |

## Examples

### Example 1: Single-residue match

**User Prompt:**
> Find "K" in "MKV".

**Expected Tool Call:**
```json
{
  "tool": "find_motif_by_pattern",
  "arguments": { "proteinSequence": "MKV", "regexPattern": "K", "motifName": "Lys" }
}
```

**Response:**
```json
{ "items": [ { "start": 1, "end": 1, "sequence": "K", "motifName": "Lys" } ] }
```

### Example 2: No match

**User Prompt:**
> Find "W" in "MKV".

**Expected Tool Call:**
```json
{
  "tool": "find_motif_by_pattern",
  "arguments": { "proteinSequence": "MKV", "regexPattern": "W" }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n) with a bounded match timeout.
- **Space Complexity:** O(number of matches).

## See Also

- [find_motif_by_prosite](find_motif_by_prosite.md) — PROSITE pattern search
- [find_protein_motifs](find_protein_motifs.md) — built-in PROSITE catalog
