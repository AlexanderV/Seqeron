# find_promoter_elements_by_matrix

Promoter elements by Bucher (1990) weight matrices (JASPAR POLII).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_promoter_elements_by_matrix` |
| **Method ID** | `MotifFinder.FindPromoterElementsByMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence with the four Bucher (1990) promoter count matrices from the JASPAR POLII collection — TATA box (POL012.1), cap signal / Inr (POL002.1), CCAAT box (POL004.1), GC box (POL003.1) — converted to log2-odds PWMs with JASPAR pseudocounts. Each matrix is scanned at the score threshold whose background false-positive rate is `falsePositiveRate` (Biopython `threshold_fpr`, precision 1000). With `bothStrands` the orientation-independent CCAAT and GC boxes are also scanned on the minus strand; TATA and cap are strand-specific. Bucher's own cut-off scores (natural-log weight scale) are not used. Order: matrix order, then position, `+` before `-`.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PromoterMatrices.cs#L81](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PromoterMatrices.cs#L81)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `falsePositiveRate` | number | No | Background false-positive rate per window and strand in [0,1] (default 0.001) |
| `bothStrands` | boolean | No | Scan CCAAT / GC boxes on both strands (default true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Hits `{name, matrixId, position, strand, sequence, score, threshold}`; position = 0-based forward window start |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | falsePositiveRate must be in [0, 1] |

## Examples

### Example 1: Synthetic promoter, FPR 1e-3

**User Prompt:**
> Find TATA / Inr / CCAAT / GC boxes in this promoter.

**Tool Call:**
```json
{
  "tool": "find_promoter_elements_by_matrix",
  "arguments": {
    "sequence": "GGGGCTATAAAAGGGGGTGGGGGCGCGTTCGTCCTCACTCTCTTCCGCATCGCTGTCTGCGAGGGCCAGCCAATCAGCGCCCCGCCCATTGGCTGGGCGGAGCC",
    "falsePositiveRate": 0.001
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "name": "TATA Box",
      "matrixId": "POL012.1",
      "position": 4,
      "strand": "+",
      "sequence": "CTATAAAAGGGGGTG",
      "score": 14.674918491724616,
      "threshold": 7.042753822501211
    },
    {
      "name": "Cap Signal",
      "matrixId": "POL002.1",
      "position": 34,
      "strand": "+",
      "sequence": "TCACTCTC",
      "score": 6.793403307766336,
      "threshold": 6.115969817666912
    },
    {
      "name": "CCAAT Box",
      "matrixId": "POL004.1",
      "position": 64,
      "strand": "+",
      "sequence": "GCCAGCCAATCA",
      "score": 12.077281742271637,
      "threshold": 7.094645943646782
    },
    {
      "name": "CCAAT Box",
      "matrixId": "POL004.1",
      "position": 85,
      "strand": "-",
      "sequence": "CCCAGCCAATGG",
      "score": 12.315014473276278,
      "threshold": 7.094645943646782
    },
    {
      "name": "GC Box",
      "matrixId": "POL003.1",
      "position": 11,
      "strand": "+",
      "sequence": "AGGGGGTGGGGGCG",
      "score": 13.369436367257704,
      "threshold": 6.836719879834334
    },
    {
      "name": "GC Box",
      "matrixId": "POL003.1",
      "position": 17,
      "strand": "+",
      "sequence": "TGGGGGCGCGTTCG",
      "score": 8.659268694378223,
      "threshold": 6.836719879834334
    },
    {
      "name": "GC Box",
      "matrixId": "POL003.1",
      "position": 76,
      "strand": "-",
      "sequence": "AATGGGCGGGGCGC",
      "score": 13.070203850439766,
      "threshold": 6.836719879834334
    }
  ]
}
```

## Worked Example

Identical hits (positions, strands, scores to 1e-5) to Biopython `pssm.search` at `threshold_fpr(1e-3)` with the same JASPAR matrices.

## See Also

- [find_regulatory_elements_both_strands](find_regulatory_elements_both_strands.md) — Consensus-pattern regulatory scan on both strands
- [pwm_score_thresholds](pwm_score_thresholds.md) — FPR thresholds for any PWM

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PromoterMatrices.cs#L81](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PromoterMatrices.cs#L81)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
