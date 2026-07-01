# find_protein_domains

Detect common protein domains via exact PROSITE patterns.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_protein_domains` |
| **Method ID** | `ProteinMotifFinder.FindDomains` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Detects common protein **domains** using EXACT PROSITE patterns: zinc finger C2H2
(PS00028 / PF00096), WD40 repeats (PS00678 / PF00400), and the protein-kinase
ATP-binding / Walker A P-loop (PS00017 / PF00069). Profile-only domains (e.g. SH3
PS50002, PDZ PS50106) have no deterministic pattern and are not detected. Each hit
reports name, Pfam accession, start, end, score and description.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L1379](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L1379)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ name, accession, start, end, score, description }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Zinc finger C2H2 (PF00096)

**User Prompt:**
> Find protein domains in "AAAACAACAAALEEEEEEEEHAAAHAAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_protein_domains",
  "arguments": { "proteinSequence": "AAAACAACAAALEEEEEEEEHAAAHAAAA" }
}
```

**Response:**
```json
{ "items": [ { "name": "Zinc Finger C2H2", "accession": "PF00096", "start": 4, "end": 24 } ] }
```
The C2H2 PROSITE consensus (PS00028) matches, from the first Cys (index 4) to the
second His (index 24).

### Example 2: No domain

**User Prompt:**
> Domains in "AAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "find_protein_domains",
  "arguments": { "proteinSequence": "AAAAAAAA" }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · D) for D catalog domains.
- **Space Complexity:** O(number of hits).

## See Also

- [find_protein_motifs](find_protein_motifs.md) — motif-level catalog scan
- [find_motif_by_prosite](find_motif_by_prosite.md) — custom PROSITE pattern
