# find_regulatory_elements_both_strands

Strand-annotated scan of the built-in regulatory element library.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_regulatory_elements_both_strands` |
| **Method ID** | `MotifFinder.FindRegulatoryElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Strand-annotated variant of `find_regulatory_elements` (`MotifFinder.FindRegulatoryElements(sequence, bothStrands: true)`): every library element (TATA, CAAT, GC box, Kozak, Shine-Dalgarno, poly(A), E-box, AP-1, NF-κB, CREB) on the given strand, plus the orientation-independent elements rescanned on the minus strand — CAAT box (Mantovani 1998), GC box (Gidoni et al. 1985) and NF-κB (enhancer, Banerji et al. 1981); AP-1, E-box and CREB are self-complementary and reported once. `position` is the 0-based forward-strand window start; `sequence` is the site read 5'→3' on its own strand. Equals Biopython `nt_search` on the sequence and its reverse complement.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.cs#L1026](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L1026)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Elements `{name, position, sequence, pattern, description, strand}` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Mixed orientation

**User Prompt:**
> Find regulatory elements on both strands.

**Tool Call:**
```json
{
  "tool": "find_regulatory_elements_both_strands",
  "arguments": {
    "sequence": "ATTGGTTTATAAACCGCCCATCCAATGGAAAGTCCCTGACTCAGGGCGGA"
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "name": "TATA Box",
      "position": 7,
      "sequence": "TATAAA",
      "pattern": "TATAAA",
      "description": "Eukaryotic core promoter element",
      "strand": "+"
    },
    {
      "name": "CAAT Box",
      "position": 0,
      "sequence": "CCAAT",
      "pattern": "CCAAT",
      "description": "Promoter element",
      "strand": "-"
    },
    {
      "name": "CAAT Box",
      "position": 21,
      "sequence": "CCAAT",
      "pattern": "CCAAT",
      "description": "Promoter element",
      "strand": "+"
    },
    {
      "name": "GC Box",
      "position": 13,
      "sequence": "GGGCGG",
      "pattern": "GGGCGG",
      "description": "Sp1 binding site",
      "strand": "-"
    },
    {
      "name": "GC Box",
      "position": 43,
      "sequence": "GGGCGG",
      "pattern": "GGGCGG",
      "description": "Sp1 binding site",
      "strand": "+"
    },
    {
      "name": "AP-1",
      "position": 36,
      "sequence": "TGACTCA",
      "pattern": "TGASTCA",
      "description": "AP-1 transcription factor binding",
      "strand": "+"
    },
    {
      "name": "NF-κB",
      "position": 26,
      "sequence": "GGGACTTTCC",
      "pattern": "GGGRNWYYCC",
      "description": "NF-κB binding site",
      "strand": "-"
    }
  ]
}
```

## See Also

- [find_regulatory_elements](find_regulatory_elements.md) — Given-strand scan
- [find_promoter_elements_by_matrix](find_promoter_elements_by_matrix.md) — Weight-matrix promoter scan

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.cs#L1026](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L1026)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
