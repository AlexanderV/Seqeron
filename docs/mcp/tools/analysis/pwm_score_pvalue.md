# pwm_score_pvalue

Exact PWM score p-value, or the exact score threshold of a p-value (Touzet & Varré 2007, TFM-Pvalue).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `pwm_score_pvalue` |
| **Method ID** | `MotifFinder.PwmScorePValue` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Exact p-value of a PWM score, P(S >= score) for a word drawn i.i.d. from the background (S = the window score exactly as scan_with_pwm / CalculatePwmScores computes it, so an observed window counts itself), or — with `pValue` — the exact score threshold: the smallest word score t with P(S >= t) <= pValue and its exact p-value. Method: Touzet & Varré (2007) TFM-Pvalue successive refinement of integer-rounded matrices (M = floor(g·W), g = 10, 100, …; error bound E), with the undecided band of words resolved by branch-and-bound enumeration of their exact scores. When the work budget is exhausted the certified bounds are returned with isExact = false and pValue = the upper bound. Give exactly one of `score` / `pValue`. Finite matrices only (use a positive pseudocount).

The inverse (`pValue`) is `MotifFinder.PwmScoreThresholdForPValue`. Unlike [pwm_score_thresholds](pwm_score_thresholds.md) (Biopython's fixed-grid approximation) no discretisation error remains: e.g. for the Wikipedia PWM, Biopython `threshold_fpr(0.01)` = 4.028388324862519 whereas the exact threshold is 4.028050165603465 (P = 0.009918212890625).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L41](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L41)
- Algorithm doc: [Position_Weight_Matrix.md](../../../algorithms/Pattern_Matching/Position_Weight_Matrix.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pwm` | object | Yes | Position Weight Matrix (finite cells) |
| `score` | number | No* | Score threshold: returns P(S >= score) |
| `pValue` | number | No* | Target p-value in [0,1]: returns the exact score threshold |
| `background` | array<number> | No | Background A,C,G,T (uniform when omitted) |

\* exactly one of `score` / `pValue`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number \| null | The requested score, or the computed threshold (null when `scoreAboveMaximum`) |
| `scoreAboveMaximum` | boolean | Even the best word has P(S >= max) > pValue: no word score qualifies, p-value 0 |
| `pValue` | number | P(S >= score); the conservative upper bound when `isExact` is false |
| `pValueLowerBound` | number | Certified lower bound |
| `pValueUpperBound` | number | Certified upper bound |
| `isExact` | boolean | The p-value is exact |
| `granularity` | number | Scale g of the integer-rounded matrix that resolved it (0 = decided without rounding) |

## Errors

| Code | Message |
|------|---------|
| 1002 | PWM is required |
| 1002 | PWM cells must be finite |
| 1003 | Background must have 4 values (A,C,G,T) |
| 1004 | Give exactly one of score or pValue |
| 1005 | Score must be finite |
| 1006 | pValue must be in [0, 1] |

## Examples

### Example 1: P-value of score 0

**User Prompt:**
> How significant is a score of 0 for this two-column PWM?

**Tool Call:**
```json
{
  "tool": "pwm_score_pvalue",
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
    "score": 0.0
  }
}
```

**Response:**
```json
{
  "score": 0,
  "scoreAboveMaximum": false,
  "pValue": 0.4375,
  "pValueLowerBound": 0.4375,
  "pValueUpperBound": 0.4375,
  "isExact": true,
  "granularity": 10
}
```

### Example 2: Threshold for p = 0.1

**Tool Call:**
```json
{
  "tool": "pwm_score_pvalue",
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
    "pValue": 0.1
  }
}
```

**Response:**
```json
{
  "score": 4,
  "scoreAboveMaximum": false,
  "pValue": 0.0625,
  "pValueLowerBound": 0.0625,
  "pValueUpperBound": 0.0625,
  "isExact": true,
  "granularity": 10
}
```

## Worked Example

Words AC score 4 (1/16), A[AGT] and [CGT]C score 0 (6/16), all others −4: P(S ≥ 0) = 7/16 = 0.4375. For p = 0.1 the smallest qualifying word score is 4 (P = 1/16); 0 would give 7/16 > 0.1.

## See Also

- [pwm_score_thresholds](pwm_score_thresholds.md) — Biopython grid-approximation thresholds
- [scan_with_pwm](scan_with_pwm.md) — Scan with the chosen threshold
- [create_pwm](create_pwm.md) — Build a PWM

## References

- Touzet H., Varré J.-S. (2007). Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2:15. Reference C++ source: CRAN `TFMPvalue` `src/Matrix.cpp`, `src/TFMpvalue.cpp`.
- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L41](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L41)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
