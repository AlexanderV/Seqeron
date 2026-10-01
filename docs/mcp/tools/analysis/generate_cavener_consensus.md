# generate_cavener_consensus

Degenerate IUPAC consensus by the Cavener (1987) rules.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_cavener_consensus` |
| **Method ID** | `MotifFinder.GenerateCavenerConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Degenerate consensus of aligned, equal-length DNA sequences by the Cavener (1987, NAR 15:1353) rules, as used by TRANSFAC and Biopython `Bio.motifs` `degenerate_consensus`. Per column, with base counts sorted c1 ≥ c2 ≥ c3 ≥ c4: a single base if c1 > 50% and c1 > 2·c2; otherwise the two-base IUPAC code if c1 + c2 > 75%; otherwise the three-base code if the fourth base is absent; otherwise `N`. Case-insensitive; only A/C/G/T are accepted (no gaps). Cross-checked against Biopython.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.cs#L491](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L491)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned DNA sequences of equal length over A/C/G/T |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | IUPAC degenerate consensus |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| — | Null element, unequal lengths or a non-ACGT character (ArgumentException) |

## Examples

### Example 1: Biopython tutorial motif

**User Prompt:**
> Cavener consensus of the Biopython tutorial instances.

**Tool Call:**
```json
{
  "tool": "generate_cavener_consensus",
  "arguments": {
    "sequences": [
      "TACAA",
      "TACGC",
      "TACAC",
      "TACCC",
      "AACCC",
      "AATGC",
      "AATGC"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "WACVC"
}
```

## Worked Example

Column 1: T×4, A×3 — T is > 50% but not > 2×A; T + A = 100% > 75% → W. Column 4: G×3, A×2, C×2, T×0 — G + A = 71% ≤ 75% and T is absent → V (A/C/G). Biopython `degenerate_consensus` returns the same WACVC.

## See Also

- [generate_consensus](generate_consensus.md) — 25%-threshold IUPAC consensus
- [generate_emboss_consensus](generate_emboss_consensus.md) — EMBOSS cons consensus

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.cs#L491](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L491)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
