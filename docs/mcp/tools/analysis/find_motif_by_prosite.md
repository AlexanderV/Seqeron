# find_motif_by_prosite

Find PROSITE-pattern matches in a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_motif_by_prosite` |
| **Method ID** | `ProteinMotifFinder.FindMotifByProsite` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts a **PROSITE pattern** to a regex internally (see `prosite_to_regex`) and then
scans a protein sequence for its overlapping matches. Each hit reports the 0-based
start/end, matched subsequence, motif name, and — as the pattern id — the original
PROSITE string, plus a heuristic score and E-value.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L265](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L265)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `prositePattern` | string | Yes | PROSITE-format pattern (min length 1) |
| `motifName` | string | No | Motif display name (default `Custom`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, sequence, motifName, pattern, score, eValue }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | PROSITE pattern cannot be null or empty |

## Examples

### Example 1: N-glycosylation-like pattern

**User Prompt:**
> Scan "MNSTV" with the PROSITE pattern "N-{P}-[ST]".

**Expected Tool Call:**
```json
{
  "tool": "find_motif_by_prosite",
  "arguments": { "proteinSequence": "MNSTV", "prositePattern": "N-{P}-[ST]", "motifName": "NGlyc" }
}
```

**Response:**
```json
{ "items": [ { "start": 1, "end": 3, "sequence": "NST", "motifName": "NGlyc", "pattern": "N-{P}-[ST]" } ] }
```
The pattern → regex `N[^P][ST]` matches NST at position 1.

### Example 2: No match

**User Prompt:**
> Scan "MKV" with "N-{P}-[ST]".

**Expected Tool Call:**
```json
{
  "tool": "find_motif_by_prosite",
  "arguments": { "proteinSequence": "MKV", "prositePattern": "N-{P}-[ST]" }
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

- [prosite_to_regex](prosite_to_regex.md) — the conversion step
- [find_motif_by_pattern](find_motif_by_pattern.md) — scan with a raw regex
