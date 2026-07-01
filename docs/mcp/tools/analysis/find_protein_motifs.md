# find_protein_motifs

Scan a protein against the built-in PROSITE-style motif catalog.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_protein_motifs` |
| **Method ID** | `ProteinMotifFinder.FindCommonMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a protein sequence against a **built-in catalog of PROSITE-style motifs**:
N-glycosylation (PS00001), PKC and CK2 phosphorylation sites, ATP/GTP P-loop, EF-hand,
zinc finger, leucine zipper, NLS, NES, SIM, WW, SH3 and more. Each hit reports its
0-based start/end, matched subsequence, motif name, pattern, a heuristic score and an
E-value.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L170](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L170)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, sequence, motifName, pattern, score, eValue }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: N-glycosylation site (PS00001)

**User Prompt:**
> Find protein motifs in "AAAANFTAAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_protein_motifs",
  "arguments": { "proteinSequence": "AAAANFTAAAA" }
}
```

**Response (relevant hit):**
```json
{ "items": [ { "start": 4, "end": 7, "sequence": "NFTA", "motifName": "ASN_GLYCOSYLATION" } ] }
```
The PS00001 pattern `N-{P}-[ST]-{P}` matches the NFTA window at index 4.

### Example 2: No N-glycosylation site

**User Prompt:**
> Protein motifs in "AAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "find_protein_motifs",
  "arguments": { "proteinSequence": "AAAAAA" }
}
```

**Response:**
```json
{ "items": [] }
```
A poly-alanine tract has no N-glycosylation (or other listed) site.

## Performance

- **Time Complexity:** O(n · P) for P catalog patterns.
- **Space Complexity:** O(number of hits).

## See Also

- [find_motif_by_prosite](find_motif_by_prosite.md) — custom PROSITE pattern search
- [find_protein_domains](find_protein_domains.md) — domain-level detection
