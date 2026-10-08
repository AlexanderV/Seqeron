# pwm_score_thresholds

PWM score thresholds (FPR / FNR / balanced / patser) from the score distribution.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `pwm_score_thresholds` |
| **Method ID** | `PositionWeightMatrix.ScoreDistribution` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds the discretised score distribution of a PWM under the motif and background models — a line-by-line port of Biopython `Bio.motifs.thresholds.ScoreDistribution` (`pssm.distribution(background, precision)`) — and returns its grid plus the derived thresholds: `thresholdFpr` (background false-positive rate ≈ `fpr`, `threshold_fpr`), `thresholdFnr` (motif false-negative rate ≈ `fnr`, `threshold_fnr`), `thresholdBalanced` with the false-positive rate it reaches (FNR ≈ FPR·`rateProportion`, `threshold_balanced`), and `thresholdPatser` (log2 FPR = −information content, `threshold_patser`). The matrix must be finite (use a positive pseudocount).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L322](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L322)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pwm` | object | Yes | Position Weight Matrix |
| `background` | array<number> | No | Background A,C,G,T (uniform when omitted) |
| `precision` | integer | No | Grid points per motif position (default 1000) |
| `fpr` | number | No | Target false-positive rate in [0,1] (default 0.01) |
| `fnr` | number | No | Target false-negative rate in [0,1] (default 0.1) |
| `rateProportion` | number | No | Balanced threshold FNR/FPR proportion (default 1.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `minScore` | number | Score of grid point 0: min(0, PWM minimum) |
| `step` | number | Grid spacing |
| `pointCount` | integer | precision × motif length |
| `meanScore` | number | Expected background score (information content) |
| `thresholdFpr` | number | Threshold for the requested FPR |
| `thresholdFnr` | number | Threshold for the requested FNR |
| `thresholdBalanced` | number | Balanced threshold |
| `balancedFalsePositiveRate` | number | FPR reached at the balanced threshold |
| `thresholdPatser` | number | patser threshold |

## Errors

| Code | Message |
|------|---------|
| 1002 | PWM is required |
| 1002 | PWM cells must be finite |
| 1003 | Background must have 4 values (A,C,G,T) |
| 1004 | Precision must be >= 1 |
| 1005 | fpr must be in [0, 1] |
| 1006 | fnr must be in [0, 1] |

## Examples

### Example 1: Two-column PWM, precision 100

**User Prompt:**
> What score threshold gives a 1% false-positive rate for this PWM?

**Tool Call:**
```json
{
  "tool": "pwm_score_thresholds",
  "arguments": {
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
    "precision": 100
  }
}
```

**Response:**
```json
{
  "minScore": -4,
  "step": 0.04020100502512563,
  "pointCount": 200,
  "meanScore": 3.25,
  "thresholdFpr": 4,
  "thresholdFnr": -0.0201005025125629,
  "thresholdBalanced": 4,
  "balancedFalsePositiveRate": 0.0625,
  "thresholdPatser": -0.0201005025125629
}
```

## Worked Example

Only the word AC scores 4 (probability 1/16 under a uniform background), so every FPR ≤ 1/16 maps to the top score 4.

## See Also

- [scan_with_pwm_both_strands](scan_with_pwm_both_strands.md) — Scan with the chosen threshold
- [create_pwm](create_pwm.md) — Build a PWM

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L322](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs#L322)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
