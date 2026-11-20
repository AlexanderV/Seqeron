# ksw2_extend

ksw2 / minimap2 extension alignment: affine or dual-affine gaps, band, Z-drop, end bonus, N penalty, CIGAR.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `ksw2_extend` |
| **Method ID** | `Ksw2Aligner` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The ksw2 / minimap2 extension-alignment functions: global or extension alignment of a query against a target with affine (q + k*e) or dual-affine (min(q + k*e, q2 + k*e2)) gap cost, a diagonal band, the Z-drop rule, the end bonus, the ambiguous-base (N) penalty and CIGAR output. variant: 'extz' (ksw_extz, scalar affine), 'extz2' (ksw_extz2_sse, minimap2 affine), 'extd' (ksw_extd, scalar dual affine), 'extd2' (ksw_extd2_sse, minimap2 -O/-E dual affine; the default). Bases are encoded A,C,G,T = 0..3 and anything else = N (4); match / mismatch / ambiguous give the matrix as minimap2 builds it (-A, -B as a penalty, -N penalty of the wildcard). Penalties are POSITIVE (minimap2 -O q,q2 -E e,e2; a gap of length k costs q + k*e). w = band half-width (-1: none; the SIMD variants need w >= 1), zdrop = off-diagonal drop-off to stop an extension (-1: none; minimap2 -z), endBonus (minimap2 -s... as the extension end bonus; SIMD variants, extensionOnly). extensionOnly = the path ends at the best cell (or at the query end when mqe + endBonus > max); otherwise it is global (corner to corner). Returned: max, zdropped, maxQ/maxT, mqe (best score reaching the query end) with mqeT, mte with mteQ, score (corner score; -1073741824 = not reached), reachEnd and the CIGAR (M, I = query only, D = target only). VERIFICATION: exact cases (no Z-drop, no approximate mode) equal an independent dynamic program on 330 cases per variant (max, mqe, mte, score; every CIGAR re-scored along its path); the Z-drop and approximate cases have not yet been compared with the real C program (it could not be run in the development environment).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/Ksw2Aligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Ksw2Aligner.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | Yes | Query sequence (ACGT, other letters = N). |
| `target` | string | Yes | Target sequence. |
| `variant` | string | No | 'extd2' (default), 'extz2', 'extd' or 'extz'. |
| `match` | integer | No | Score of a match (minimap2 -A). (default: 2) |
| `mismatch` | integer | No | Penalty of a mismatch, positive (minimap2 -B). (default: 4) |
| `ambiguous` | integer | No | Penalty of the ambiguous base N, positive (minimap2 -N / sc_ambi); 0 gives -e. (default: 1) |
| `gapOpen` | integer | No | Gap open penalty q (minimap2 -O first value), >= 0. (default: 4) |
| `gapExtend` | integer | No | Gap extension penalty e (minimap2 -E first value), >= 1. (default: 2) |
| `gapOpen2` | integer | No | Second gap open penalty q2 (dual variants; minimap2 -O second value). (default: 24) |
| `gapExtend2` | integer | No | Second gap extension penalty e2 (dual variants; minimap2 -E second value). (default: 1) |
| `band` | integer | No | Band half-width (-1 = none). (default: -1) |
| `zdrop` | integer | No | Z-drop (-1 = none; minimap2 -z). (default: -1) |
| `endBonus` | integer | No | End bonus (SIMD variants, extensionOnly). (default: 0) |
| `extensionOnly` | boolean | No | Extension alignment: the path ends at the best cell (or the query end with the end bonus). (default: false) |
| `rightAlignGaps` | boolean | No | Right-align gaps (KSW_EZ_RIGHT). (default: false) |
| `scoreOnly` | boolean | No | Score only, no CIGAR. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `max` | integer | Best score of the extension (0 when nothing scored). |
| `zdropped` | boolean | The extension stopped by the Z-drop rule or an empty band. |
| `maxQ` | integer | Query coordinate of the best cell (-1 none). |
| `maxT` | integer | Target coordinate of the best cell (-1 none). |
| `mqe` | integer | Best score reaching the query end. |
| `mqeT` | integer | Target coordinate of mqe. |
| `mte` | integer | Best score reaching the target end. |
| `mteQ` | integer | Query coordinate of mte. |
| `score` | integer | Corner score (-1073741824 = not reached). |
| `reachEnd` | boolean | The path was taken to the query end because mqe + endBonus > max. |
| `cigar` | string | CIGAR (M, I = query only, D = target only; * when empty). |

## Errors

| Code | Message |
|------|---------|
| 1001 | query / target is null (an empty sequence aligns to nothing: max 0, cigar *) |
| 1002 | variant unknown; gapOpen < 0 or gapExtend < 1 (byte range) |

## Examples

### Example 1: Identical 8-mers

**Tool Call:**
```json
{
  "tool": "ksw2_extend",
  "arguments": {
    "query": "ACGTACGT",
    "target": "ACGTACGT"
  }
}
```

**Response:**
```json
{
  "max": 16,
  "zdropped": false,
  "maxQ": 7,
  "maxT": 7,
  "mqe": 16,
  "mqeT": 7,
  "mte": 16,
  "mteQ": 7,
  "score": 16,
  "reachEnd": false,
  "cigar": "8M"
}
```

## See Also

- [write_alignment](write_alignment.md)
- [stretcher_align](stretcher_align.md)
