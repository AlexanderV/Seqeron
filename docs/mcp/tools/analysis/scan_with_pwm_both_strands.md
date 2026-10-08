# scan_with_pwm_both_strands

Scan both strands of a DNA sequence with a PWM (Biopython search(both=True)).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `scan_with_pwm_both_strands` |
| **Method ID** | `MotifFinder.ScanWithPwmBothStrands` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans both strands with a 4×L Position Weight Matrix (rows A,C,G,T) as Biopython `PositionSpecificScoringMatrix.search(sequence, threshold, both=True)`: the plus strand is scored with the PWM and the minus strand with its reverse complement over the same forward windows; every hit with score ≥ `threshold` is reported. `position` is the 0-based forward-strand window start on both strands; `biopythonPosition` is Biopython's coordinate (minus strand: start − n). The minus-strand `matchedSequence` is the site read 5'→3' on the minus strand. Order: ascending position, `+` before `-`. Windows with non-ACGT symbols are skipped.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L100](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L100)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan |
| `pwm` | object | Yes | Position Weight Matrix |
| `threshold` | number | No | Minimum score (inclusive, default 0.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Hits `{position, biopythonPosition, strand, matchedSequence, pattern, score}` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | PWM is required |
| 1002 | Matrix must have 4 rows |

## Examples

### Example 1: Two-column PWM for AC

**User Prompt:**
> Scan ACGTAC on both strands with a PWM favouring AC.

**Tool Call:**
```json
{
  "tool": "scan_with_pwm_both_strands",
  "arguments": {
    "sequence": "ACGTAC",
    "pwm": {
      "matrix": [
        [
          2.0,
          -2.0
        ],
        [
          -2.0,
          2.0
        ],
        [
          -2.0,
          -2.0
        ],
        [
          -2.0,
          -2.0
        ]
      ],
      "length": 2
    },
    "threshold": 0.0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "position": 0,
      "biopythonPosition": 0,
      "strand": "+",
      "matchedSequence": "AC",
      "pattern": "AC",
      "score": 4
    },
    {
      "position": 2,
      "biopythonPosition": -4,
      "strand": "-",
      "matchedSequence": "AC",
      "pattern": "AC",
      "score": 4
    },
    {
      "position": 4,
      "biopythonPosition": 4,
      "strand": "+",
      "matchedSequence": "AC",
      "pattern": "AC",
      "score": 4
    }
  ]
}
```

## Worked Example

The forward window GT at position 2 reads AC on the minus strand, so it scores 4 there (Biopython position 2 − 6 = −4).

## See Also

- [scan_with_pwm](scan_with_pwm.md) — Forward-strand scan
- [pwm_score_thresholds](pwm_score_thresholds.md) — Choose a threshold from the score distribution
- [create_pwm](create_pwm.md) — Build a PWM

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L100](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L100)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
