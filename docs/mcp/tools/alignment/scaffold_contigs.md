# scaffold_contigs

Join contigs into scaffolds using paired-end link records, inserting a run of gapCharacter (length = gapSize; unknown/non-positive -> 100) between linked contigs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `scaffold_contigs` |
| **Method ID** | `SequenceAssembler.Scaffold` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Join contigs into scaffolds using paired-end link records, inserting a run of gapCharacter (length = gapSize; unknown/non-positive -> 100) between linked contigs.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L679](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L679)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contigs` | array<string> | Yes | Contig sequences. |
| `links` | array<object> | Yes | Paired-end links between contigs, each an object {contig1, contig2, gapSize}. |
| `gapCharacter` | string | No | Single-character gap filler (e.g. "N"). (default: N) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `scaffolds` | array<string> | The assembled scaffolds, one string per scaffold. |

## Errors

| Code | Message |
|------|---------|
| 1013 | gapCharacter must be exactly one character |

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

## Worked Example

Link 0 -> 1 with gapSize 3 places contig[1] after contig[0] separated by three `N`s: `AAA` + `NNN` + `GGG`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L679](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L679)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
