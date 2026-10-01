# generate_emboss_consensus

EMBOSS cons scoring-matrix plurality consensus of an alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_emboss_consensus` |
| **Method ID** | `MotifFinder.GenerateEmbossConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Line-by-line port of EMBOSS 6.6.0 `cons` (`embConsCalc`): per column the residue with the highest positive-match score under the scoring matrix (EDNAFULL for nucleotides, EBLOSUM62 for proteins) is chosen; if its positive-match weight is below `plurality` (default half the total weight) or fewer than `identity` residues are identical, the no-consensus symbol (`N` nucleotide / `X` protein) is emitted; residues whose weight is ≤ `setcase` are written in lower case. Gaps `-`, `.`, `~` are allowed. Single-precision arithmetic as in EMBOSS; outputs are character-identical to `cons -plurality P -identity I -setcase S` (780/780 reference runs). `residueType: "auto"` types the alignment exactly as `cons` does without `-snucleotide`/`-sprotein` (each sequence typed on reading, the set takes the first sequence's type — matrix and N/X; 3,700/3,700 reference runs), and `padRaggedRows: true` pads shorter rows with trailing gaps as `cons` does (`ajSeqsetFill`). With the explicit `protein` type the no-consensus symbol is always `X` (`cons -sprotein` takes N/X from the first sequence).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L152](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L152)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least two aligned sequences (DNA or protein); equal length unless `padRaggedRows` |
| `residueType` | string | No | 'nucleotide' (EDNAFULL, default), 'protein' (EBLOSUM62) or 'auto' (decided from the first sequence like cons) |
| `plurality` | number | No | Minimum positive-match weight (cons -plurality); default half the total sequence weight |
| `identity` | integer | No | Required number of identical residues (cons -identity, default 0 = off) |
| `setcase` | number | No | Weight at or below which output is lower case (cons -setcase); default half the total sequence weight |
| `weights` | array<number> | No | Optional per-sequence weights (finite, >= 0; default 1.0 each) |
| `padRaggedRows` | boolean | No | Pad rows shorter than the longest with trailing '-' (cons ajSeqsetFill) instead of rejecting them (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | EMBOSS cons consensus (lower case where weight ≤ setcase) |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least two sequences are required |
| 1002 | residueType must be 'nucleotide', 'protein' or 'auto' |
| — | Null element, unequal lengths, invalid character or weight count (ArgumentException); negative identity / invalid weight (ArgumentOutOfRangeException) |

## Examples

### Example 1: cons -snucleotide

**User Prompt:**
> EMBOSS cons consensus of this DNA alignment.

**Tool Call:**
```json
{
  "tool": "generate_emboss_consensus",
  "arguments": {
    "sequences": [
      "ACGTAC-T",
      "ACGTTCAT",
      "AGGTAC-T",
      "tCGAAG-T"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "ACGTACnT"
}
```

### Example 2: cons -sprotein -plurality 1.0 -setcase 3.0

**Tool Call:**
```json
{
  "tool": "generate_emboss_consensus",
  "arguments": {
    "sequences": [
      "MKVLAAGIVG",
      "MKVLSAGIVA",
      "MRVLAAG-VG",
      "MKILTAGLVG"
    ],
    "residueType": "protein",
    "plurality": 1.0,
    "setcase": 3.0
  }
}
```

**Response:**
```json
{
  "consensus": "MKVLaAGiVg"
}
```

### Example 3: cons without a type flag, ragged rows

**Tool Call:**
```json
{
  "tool": "generate_emboss_consensus",
  "arguments": {
    "sequences": ["ACGTAC", "ACG", "AC"],
    "residueType": "auto",
    "padRaggedRows": true
  }
}
```

**Response:**
```json
{
  "consensus": "ACGnnn"
}
```

## See Also

- [generate_dumb_consensus](generate_dumb_consensus.md) — Biopython dumb_consensus
- [generate_consensus](generate_consensus.md) — 25%-threshold IUPAC consensus

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L152](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L152)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
