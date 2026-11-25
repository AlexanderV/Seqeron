# find_acceptor_sites

Find candidate 3' (acceptor) splice sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_acceptor_sites` |
| **Method ID** | `SpliceSitePredictor.FindAcceptorSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence (upper-cased, T→U) for 3' splice acceptor sites at every AG dinucleotide.

- **`model = "shapiro-senapathy"` (default):** the Shapiro & Senapathy (1987) 3' ss score of the 15-nt
  subsequence -14..+1 with the Table 1 subtable `category`:
  `((t1 − l1)/(h1 − l1) + (t2 − l2)/(h2 − l2)) / 2`, where the pyrimidine region -14..-5 contributes the
  best 8 of its 10 percentages and the CAGG part the percentages at -3, -2, -1, +1 (position -4 unscored);
  i.e. the published 0-100 score / 100 (reproduces the paper's Table 5 3' ss scores, 40/41 at 1 dp). A
  candidate with an AG inside -14..-5 is eliminated (as in the paper); a window that runs off the sequence
  or contains a non-ACGTU character is not scored. Confidence = clamp((score − 0.5)/0.5, 0, 1).
- **`model = "maxent"`:** MaxEntScan score3 (Yeo & Burge 2004) of the 23-nt window -20..+3, in bits
  (`minScore` in bits); confidence = 2^s/(1+2^s).
- **`includeNonCanonical` (shapiro-senapathy only):** U12-type acceptors (AC or AG terminus), scored as the
  best match fraction to the U12 branch site `UUCCUUAAC` (branch A = second A; Brock et al. 2008) with the
  branch A 8-21 nt upstream of the last intron nt (Levine & Durbin 2001); type `U12Acceptor`.

`position` is the 0-based index of the G of AG (last intron nucleotide); `motif` is the -20..+3 context
(truncated at the sequence ends).

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L579](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L579)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `minScore` | number | No | Minimum score (default 0.5): [0,1] for shapiro-senapathy / U12, bits for maxent |
| `includeNonCanonical` | boolean | No | Also report U12-type acceptors (default false) |
| `model` | string | No | `shapiro-senapathy` (default) or `maxent` |
| `category` | string | No | S&S Table 1 subtable: `primate` (default), `total`, `rodent`, `othermammal`, `othervertebrate`, `invertebrate`, `plant`, `viral` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based index of the G of the AG (last intron nt) |
| `sites[].type` | string | `Acceptor` or `U12Acceptor` |
| `sites[].motif` | string | -20..+3 context (RNA alphabet) |
| `sites[].score` | number | S&S score [0,1], MaxEntScan bits, or U12 branch-site match fraction |
| `sites[].confidence` | number | Confidence [0,1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown model '...'; expected 'shapiro-senapathy' or 'maxent'. |
| 1003 | Unknown category '...' |

## Examples

### Example 1: Published primate 3' ss (S&S 1987 Table 5, HLA intron 6, 95.8)

`TTATTCTACTCCAGG`:

**Response:**
```json
{ "sites": [ { "position": 13, "type": "Acceptor", "motif": "UUAUUCUACUCCAGG", "score": 0.9578125, "confidence": 0.915625 } ] }
```

### Example 2: MaxEntScan model

`TTCCAAACGAACTTTTGTAGGGA`, `model = "maxent"`, `minScore = 0` → position 19, score 2.8867730651152104 bits
(MaxEntScan score3 reference 2.89).

### Example 3: No AG dinucleotide

```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_donor_sites](find_donor_sites.md) — 5' donor sites
- [find_branch_points](find_branch_points.md) — intronic branch points
- [predict_introns](predict_introns.md) — pair donors and acceptors into introns
