# alphabet_pwm_score_pvalue

Exact p-value of a score, or the exact score threshold of a p-value, for a PWM over any alphabet (protein, RNA, …).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `alphabet_pwm_score_pvalue` |
| **Method ID** | `MotifFinder.AlphabetPwmScorePValue` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds a log-odds PWM over an arbitrary alphabet from aligned instances (as [create_alphabet_pwm](create_alphabet_pwm.md), Biopython `motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)`; the pseudocount must make every cell finite, default 0.5) and returns the exact p-value P(S >= score) for a word drawn i.i.d. from the background (S = the window score of [scan_with_alphabet_pwm](scan_with_alphabet_pwm.md)), or — with `pValue` — the smallest word score t with P(S >= t) <= pValue and its exact p-value. Same engine as [pwm_score_pvalue](pwm_score_pvalue.md) (Touzet & Varré 2007 TFM-Pvalue refinement + exact band enumeration) with K = |alphabet| rows; same search options. The inverse is `MotifFinder.AlphabetPwmScoreThresholdForPValue`. Give exactly one of `score` / `pValue`.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L190](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L190)
- Algorithm doc: [Position_Weight_Matrix.md](../../../algorithms/Pattern_Matching/Position_Weight_Matrix.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned instances of equal length |
| `alphabet` | string | Yes | Alphabet in row order (e.g. ACDEFGHIKLMNPQRSTVWY) |
| `score` | number | No* | Score threshold: returns P(S >= score) |
| `pValue` | number | No* | Target p-value in [0,1]: returns the exact score threshold |
| `pseudocount` | number | No | Pseudocount per cell (default 0.5) |
| `pseudocounts` | array<number> | No | Per-symbol pseudocounts in alphabet order |
| `background` | array<number> | No | Background in alphabet order (log-odds and word distribution; uniform when omitted) |
| `ignoreUnknownSymbols` | boolean | No | Skip symbols outside the alphabet when counting (default false) |
| `initialGranularity` | number | No | TFM-Pvalue initial granularity (default 0.1) |
| `maxGranularity` | number | No | TFM-Pvalue finest granularity |
| `decreaseFactor` | number | No | Granularity decrease factor (default 10) |
| `maxStates` | integer | No | DP states per column (default 2097152) |
| `maxSuffixSet` | integer | No | Suffix-score set size (default 1048576) |
| `exhaustive` | boolean | No | Always exact (default false) |

\* exactly one of `score` / `pValue`.

## Output Schema

Same as [pwm_score_pvalue](pwm_score_pvalue.md): `score` (null when `scoreAboveMaximum`), `scoreAboveMaximum`, `pValue`, `pValueLowerBound`, `pValueUpperBound`, `isExact`, `granularity`.

## Errors

| Code | Message |
|------|---------|
| 1002 | At least one sequence is required. |
| 1002 | Alphabet is required (e.g. ACDEFGHIKLMNPQRSTVWY). |
| 1003 | PWM cells must be finite: use a positive pseudocount. |
| 1003 | background must have one value per alphabet symbol. |
| 1004 | Give exactly one of score or pValue |
| 1005 | Score must be finite |
| 1006 | pValue must be in [0, 1] |

## Examples

### Example 1: P-value of a protein site score

**User Prompt:**
> How unusual is a score of 12.94 (the MRVLGT instance) for this protein motif?

**Tool Call:**
```json
{
  "tool": "alphabet_pwm_score_pvalue",
  "arguments": {
    "sequences": ["MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT"],
    "alphabet": "ACDEFGHIKLMNPQRSTVWY",
    "score": 12.939220045315675
  }
}
```

**Response:**
```json
{
  "score": 12.939220045315675,
  "scoreAboveMaximum": false,
  "pValue": 1.953125000000001E-06,
  "pValueLowerBound": 1.953125000000001E-06,
  "pValueUpperBound": 1.953125000000001E-06,
  "isExact": true,
  "granularity": 10
}
```

### Example 2: Threshold for p = 1e-5

**Tool Call:**
```json
{
  "tool": "alphabet_pwm_score_pvalue",
  "arguments": {
    "sequences": ["MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT"],
    "alphabet": "ACDEFGHIKLMNPQRSTVWY",
    "pValue": 1e-5
  }
}
```

**Response:**
```json
{
  "score": 11.064750927399535,
  "scoreAboveMaximum": false,
  "pValue": 8.906250000000003E-06,
  "pValueLowerBound": 8.906250000000003E-06,
  "pValueUpperBound": 8.906250000000003E-06,
  "isExact": true,
  "granularity": 10
}
```

## Worked Example

Exhaustive enumeration of all 20^6 = 64,000,000 hexapeptides (uniform background): P(S ≥ 12.939220045315675) = 1.953125e-06; the smallest word score with P ≤ 1e-5 is 11.064750927399535 (P = 8.90625e-06); the next lower word score 11.064750927399533 has P = 1.21875e-05.

## See Also

- [alphabet_pwm_score_thresholds](alphabet_pwm_score_thresholds.md) — Biopython grid-approximation thresholds
- [scan_with_alphabet_pwm](scan_with_alphabet_pwm.md) — Scan with the chosen threshold
- [pwm_score_pvalue](pwm_score_pvalue.md) — DNA PWMs, Markov backgrounds

## References

- Touzet H., Varré J.-S. (2007). Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2:15.
- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L190](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L190)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
