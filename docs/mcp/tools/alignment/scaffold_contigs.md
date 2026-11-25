# scaffold_contigs

Order and orient contigs into scaffolds from paired-end links with the greedy path-merging rule (Huson, Reinert & Myers 2002): links in descending weight, a link used only if both contig ends it joins are free and it closes no cycle. Positive gapSize -> that many gapCharacter; negative gapSize -> the -gapSize overlap is collapsed when its bases agree (N matches anything), else (and for 0) the AGP unknown gap of 100.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `scaffold_contigs` |
| **Method ID** | `SequenceAssembler.Scaffold` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Order and orient contigs into scaffolds from paired-end links with the greedy path-merging rule (Huson, Reinert & Myers 2002): links in descending weight, a link used only if both contig ends it joins are free and it closes no cycle. Positive gapSize -> that many gapCharacter; negative gapSize -> the -gapSize overlap is collapsed when its bases agree (N matches anything), else (and for 0) the AGP unknown gap of 100.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1121)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contigs` | array<string> | Yes | Contig sequences. |
| `links` | array<object> | Yes | Links `{contig1, contig2, gapSize, contig1Reversed?, contig2Reversed?, weight?}`: contig2 (reverse-complemented if contig2Reversed) follows contig1 (reverse-complemented if contig1Reversed) after gapSize characters (negative = overlap); weight = supporting mate pairs, >= 1 (default 1). |
| `gapCharacter` | string | No | Single-character gap filler (e.g. "N"). (default: N) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `scaffolds` | array<string> | The assembled scaffolds, one string per scaffold. |

## Errors

| Code | Message |
|------|---------|
| 1013 | gapCharacter must be exactly one character |
| — | ArgumentOutOfRangeException: a link weight < 1 |

## Examples

### Example 1: One link

**Tool Call:**
```json
{
  "tool": "scaffold_contigs",
  "arguments": {
    "contigs": [
      "AAA",
      "GGG"
    ],
    "links": [
      {
        "contig1": 0,
        "contig2": 1,
        "gapSize": 3
      }
    ],
    "gapCharacter": "N"
  }
}
```

**Response:**
```json
{
  "scaffolds": [
    "AAANNNGGG"
  ]
}
```

### Example 2: No links

**Tool Call:**
```json
{
  "tool": "scaffold_contigs",
  "arguments": {
    "contigs": [
      "AAA",
      "GGG"
    ],
    "links": [],
    "gapCharacter": "N"
  }
}
```

**Response:**
```json
{
  "scaffolds": [
    "AAA",
    "GGG"
  ]
}
```

### Example 3: Orientation, weight and overlap

**Tool Call:**
```json
{
  "tool": "scaffold_contigs",
  "arguments": {
    "contigs": ["AAAA", "CCCC", "GGTT"],
    "links": [
      { "contig1": 0, "contig2": 1, "gapSize": 2 },
      { "contig1": 0, "contig2": 2, "gapSize": 3, "contig2Reversed": true, "weight": 4 }
    ]
  }
}
```

**Response:**
```json
{
  "scaffolds": ["AAAANNNAACC", "CCCC"]
}
```

The weight-4 link is taken first and joins the 3' end of contig 0 to contig 2 read reverse-complemented (`AACC`); the 0 -> 1 link then needs the same (already joined) end and is rejected, so contig 1 is a singleton. With `contigs = ["TACGGA", "ACGTAC"]` and a link `{1, 0, -3}`, the 3-base overlap `TAC` agrees and is collapsed: `ACGTACGGA`.

## Worked Example

Link 0 -> 1 with gapSize 3 places contig[1] after contig[0] separated by three `N`s: `AAA` + `NNN` + `GGG`. Algorithm and sources: [Scaffolding.md](../../../algorithms/Extended_Assembly/Scaffolding.md).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1121)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
