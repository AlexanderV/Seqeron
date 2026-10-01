# scan_with_alphabet_pwm

Score and scan a sequence with a PWM over any alphabet.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `scan_with_alphabet_pwm` |
| **Method ID** | `MotifFinder.ScanWithAlphabetPwm` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds the PWM of `create_alphabet_pwm` from the aligned instances, then returns the score of every window (generic form of Biopython `pssm.calculate`: Σ_j W[s[i+j], j]; a window containing a symbol outside the alphabet scores NaN, the `_pwm.c` rule — Biopython 1.88 `calculate` itself only accepts DNA) and the forward-strand hits with score ≥ `threshold` in ascending position (Biopython `search(sequence, threshold, both=False)`). Case-insensitive. `scores` holds null where the score is not finite; `invalidWindows` lists the NaN windows (any other null is −∞).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L142](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L142)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to scan |
| `sequences` | array<string> | Yes | Aligned instances defining the motif |
| `alphabet` | string | Yes | Alphabet in row order |
| `threshold` | number | No | Minimum hit score, finite (default 0) |
| `pseudocount` | number | No | Pseudocount per cell (default 0) |
| `pseudocounts` | array<number> | No | Per-symbol pseudocounts in alphabet order |
| `background` | array<number> | No | Background in alphabet order (uniform when omitted) |
| `ignoreUnknownSymbols` | boolean | No | Skip symbols outside the alphabet when counting the instances (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | PWM consensus |
| `length` | integer | Motif length |
| `scores` | array<number|null> | Score of every window i = 0 … n − L; null when not finite |
| `invalidWindows` | array<integer> | Windows containing a symbol outside the alphabet (NaN) |
| `hits` | array<object> | Hits `{position, matchedSequence, pattern, score}` |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required. |
| 1001 | Alphabet is required (e.g. ACDEFGHIKLMNPQRSTVWY). |
| 1002 | threshold must be finite. |

## Examples

### Example 1: Protein motif with unknown symbols

**User Prompt:**
> Scan GGMKVLATxxMRvLGTPPMKXLAS with the motif MKVLAT/MRVLGT/MKILAS/LKVLAT/MKVIAT/MEVLAT (protein alphabet, pseudocount 0.5), hits ≥ 5.

**Tool Call:**
```json
{
  "tool": "scan_with_alphabet_pwm",
  "arguments": {
    "sequence": "GGMKVLATxxMRvLGTPPMKXLAS",
    "sequences": [
      "MKVLAT",
      "MRVLGT",
      "MKILAS",
      "LKVLAT",
      "MKVIAT",
      "MEVLAT"
    ],
    "alphabet": "ACDEFGHIKLMNPQRSTVWY",
    "threshold": 5.0,
    "pseudocount": 0.5
  }
}
```

**Response:**
```json
{
  "consensus": "MKVLAT",
  "length": 6,
  "scores": [
    -4.068431430675826,
    -4.068431430675826,
    16.398651663952972,
    null,
    null,
    null,
    null,
    null,
    null,
    null,
    12.939220045315675,
    -4.068431430675826,
    -4.068431430675826,
    -2.4834689299546704,
    -4.068431430675826,
    null,
    null,
    null,
    null
  ],
  "invalidWindows": [
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    15,
    16,
    17,
    18
  ],
  "hits": [
    {
      "position": 2,
      "matchedSequence": "MKVLAT",
      "pattern": "MKVLAT",
      "score": 16.398651663952972
    },
    {
      "position": 10,
      "matchedSequence": "MRvLGT",
      "pattern": "MKVLAT",
      "score": 12.939220045315675
    }
  ]
}
```

## Worked Example

Window scores equal Σ_j `pssm[letter][j]` of the Biopython 1.88 PSSM (≤ 1e-14); the lower-case `v` is scored as V; windows touching `x` / `X` are NaN. For alphabet ACGT the scores equal Biopython `calculate` (float32) exactly and the hits equal `search(both=False)` (92 random DNA cases).

## See Also

- [create_alphabet_pwm](create_alphabet_pwm.md) — Build the matrix and its statistics
- [scan_with_pwm](scan_with_pwm.md) — DNA PWM scan

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L142](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L142)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
