# alphabet_pwm_score_thresholds

Score thresholds of a PWM over any alphabet from its discretised score distribution (Biopython `pssm.distribution`).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `alphabet_pwm_score_thresholds` |
| **Method ID** | `AlphabetPositionWeightMatrix.ScoreDistribution` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds a log-odds PWM over an arbitrary alphabet from aligned instances (as [create_alphabet_pwm](create_alphabet_pwm.md); positive pseudocount, default 0.5) and returns the thresholds of Biopython's `pssm.distribution(background, precision)` (`Bio.motifs.thresholds.ScoreDistribution`, N. Dojer 2008): the grid (min score, step, points = precision · L), the background FPR threshold, motif FNR threshold, balanced threshold (FNR = FPR · rateProportion) and patser threshold (log2 FPR = −information content). Biopython's DP iterates the PSSM's own alphabet, so it applies unchanged to protein PSSMs; this is the K-row form of [pwm_score_thresholds](pwm_score_thresholds.md). For exact p-values use [alphabet_pwm_score_pvalue](alphabet_pwm_score_pvalue.md).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L394](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L394)
- Algorithm doc: [Position_Weight_Matrix.md](../../../algorithms/Pattern_Matching/Position_Weight_Matrix.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned instances of equal length |
| `alphabet` | string | Yes | Alphabet in row order |
| `pseudocount` | number | No | Pseudocount per cell (default 0.5) |
| `pseudocounts` | array<number> | No | Per-symbol pseudocounts in alphabet order |
| `background` | array<number> | No | Background in alphabet order (uniform when omitted) |
| `ignoreUnknownSymbols` | boolean | No | Skip symbols outside the alphabet when counting (default false) |
| `precision` | integer | No | Grid points per motif position (default 1000) |
| `fpr` | number | No | Target background false-positive rate (default 0.01) |
| `fnr` | number | No | Target motif false-negative rate (default 0.1) |
| `rateProportion` | number | No | Balanced threshold FNR/FPR (default 1.0) |

## Output Schema

Same as [pwm_score_thresholds](pwm_score_thresholds.md): `minScore`, `step`, `pointCount`, `meanScore`, `thresholdFpr`, `thresholdFnr`, `thresholdBalanced`, `balancedFalsePositiveRate`, `thresholdPatser`.

## Errors

| Code | Message |
|------|---------|
| 1002 | At least one sequence is required. |
| 1003 | PWM cells must be finite: use a positive pseudocount. |
| 1003 | Precision must be >= 1. |
| 1004 | fpr must be in [0, 1]. |
| 1004 | fnr must be in [0, 1]. |

## Examples

### Example 1: Protein motif thresholds (Biopython defaults)

**Tool Call:**
```json
{
  "tool": "alphabet_pwm_score_thresholds",
  "arguments": {
    "sequences": ["MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT"],
    "alphabet": "ACDEFGHIKLMNPQRSTVWY"
  }
}
```

**Response:**
```json
{
  "minScore": -4.068431430675826,
  "step": 0.0034117491406282377,
  "pointCount": 6000,
  "meanScore": 3.809139711610947,
  "thresholdFpr": 2.843772328236984,
  "thresholdFnr": -0.6157413003600496,
  "thresholdBalanced": -0.6157413003600496,
  "balancedFalsePositiveRate": 0.23819500000000032,
  "thresholdPatser": 0.970722050032081
}
```

## Worked Example

Biopython 1.88: `motifs.create(instances, alphabet="ACDEFGHIKLMNPQRSTVWY").counts.normalize(pseudocounts=0.5).log_odds().distribution(precision=1000)` gives `threshold_fpr(0.01)` = 2.843772328236984, `threshold_fnr(0.1)` = −0.6157413003600496, `threshold_balanced()` = −0.6157413003600496, `threshold_patser()` = 0.970722050032081 — identical.

## See Also

- [alphabet_pwm_score_pvalue](alphabet_pwm_score_pvalue.md) — Exact p-values and thresholds
- [pwm_score_thresholds](pwm_score_thresholds.md) — DNA PWMs
- [create_alphabet_pwm](create_alphabet_pwm.md) — Build the PWM

## References

- Biopython `Bio.motifs.thresholds.ScoreDistribution` (N. Dojer 2008, adapted by B. Wilczynski); `PositionSpecificScoringMatrix.distribution`.
- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L394](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L394)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
