# find_sigma70_promoters

Pair bacterial σ70 −35 / −10 consensus boxes over a 15–21 bp spacer.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_sigma70_promoters` |
| **Method ID** | `MotifFinder.FindSigma70Promoters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairs σ70 promoter boxes per Harley & Reynolds (1987, NAR 15:2343): every −35 hexamer within `maxMismatches35` of TTGACA followed, after a spacer of `minSpacer` … `maxSpacer` bp (15–21; 17 ± 1 in 92 % of promoters), by a −10 hexamer within `maxMismatches10` of TATAAT. Each candidate reports both mismatch counts, their sum and |spacer − 17|. Coordinates are 0-based forward-strand starts; with `bothStrands` the reverse complement is scanned too (minus-strand boxes read 5'→3' on that strand). Order: `+` before `-`, then −35 start, then spacer. For a quantitative promoter score use `predict_sigma70_promoters`.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L47](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L47)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `maxMismatches35` | integer | No | Max mismatches to TTGACA, 0–6 (default 2) |
| `maxMismatches10` | integer | No | Max mismatches to TATAAT, 0–6 (default 2) |
| `minSpacer` | integer | No | Minimum spacer (default 15) |
| `maxSpacer` | integer | No | Maximum spacer (default 21) |
| `bothStrands` | boolean | No | Also scan the minus strand (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Candidates `{strand, minus35Start, minus35, minus10Start, minus10, spacer, mismatches35, mismatches10, totalMismatches, spacerDeviation}` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Must be in 0 … 6. |
| 1002 | Must be ≥ minSpacer. |

## Examples

### Example 1: E. coli lacUV5 promoter

**User Prompt:**
> Find σ70 promoters in the lac control region (≤ 1 mismatch in −35, exact −10).

**Tool Call:**
```json
{
  "tool": "find_sigma70_promoters",
  "arguments": {
    "sequence": "TCAGCATTCGAGCTTACGGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT",
    "maxMismatches35": 1,
    "maxMismatches10": 0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "strand": "+",
      "minus35Start": 68,
      "minus35": "TTTACA",
      "minus10Start": 92,
      "minus10": "TATAAT",
      "spacer": 18,
      "mismatches35": 1,
      "mismatches10": 0,
      "totalMismatches": 1,
      "spacerDeviation": 1
    }
  ]
}
```

## Worked Example

The lacUV5 −35 TTTACA (1 mismatch) at 68 pairs with the −10 TATAAT at 92 over an 18-bp spacer. 300 random sequences with planted boxes (358 candidates, both strands, random limits) are identical to an independent Python brute force.

## See Also

- [predict_sigma70_promoters](predict_sigma70_promoters.md) — Promoter Calculator free-energy model
- [find_regulatory_elements](find_regulatory_elements.md) — Consensus regulatory-element library

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L47](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L47)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
