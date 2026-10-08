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

Exact p-value of a PWM score, P(S >= score) for a word drawn from the background — i.i.d. (`background`) or an order-m Markov chain (`markovFrequencies`, an RSAT (m+1)-mer frequency table: P(w) = P(w[0..m−1])·∏P(w_c | previous m letters), `MotifFinder.PwmMarkovScorePValue`; the DP state carries the last m letters, as MACRO-APE's dinucleotide-background distribution) (S = the window score exactly as scan_with_pwm / CalculatePwmScores computes it, so an observed window counts itself), or — with `pValue` — the exact score threshold: the smallest word score t with P(S >= t) <= pValue and its exact p-value. Method: Touzet & Varré (2007) TFM-Pvalue successive refinement of integer-rounded matrices (M = floor(g·W), g = 10, 100, …; error bound E), with the undecided band of words resolved by branch-and-bound enumeration of their exact scores. When the work budget is exhausted the certified bounds are returned with isExact = false and pValue = the upper bound. Give exactly one of `score` / `pValue`. Finite matrices only (use a positive pseudocount).

Search options mirror the TFM-Pvalue driver (`initialGranularity` 0.1, `maxGranularity`, `decreaseFactor` 10 — granularity = rounding step 1/g; the reported `granularity` is the scale g) plus the work budgets `maxStates` (2^21) / `maxSuffixSet` (2^20); `exhaustive: true` resolves any band left undecided by unbounded enumeration, so the result is always exact (time/memory permitting). Defaults reproduce the plain call bit for bit.

The inverse (`pValue`) is `MotifFinder.PwmScoreThresholdForPValue`. Unlike [pwm_score_thresholds](pwm_score_thresholds.md) (Biopython's fixed-grid approximation) no discretisation error remains: e.g. for the Wikipedia PWM, Biopython `threshold_fpr(0.01)` = 4.028388324862519 whereas the exact threshold is 4.028050165603465 (P = 0.009918212890625).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L45](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L45)
- Algorithm doc: [Position_Weight_Matrix.md](../../../algorithms/Pattern_Matching/Position_Weight_Matrix.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pwm` | object | Yes | Position Weight Matrix (finite cells) |
| `score` | number | No* | Score threshold: returns P(S >= score) |
| `pValue` | number | No* | Target p-value in [0,1]: returns the exact score threshold |
| `background` | array<number> | No | i.i.d. background A,C,G,T (uniform when omitted); not with `markovFrequencies` |
| `markovFrequencies` | object | No | Markov background: (m+1)-mer → frequency (RSAT `-bgfile` oligos; m ≤ 10) |
| `markovPseudoFrequency` | number | No | RSAT pseudo-frequency ψ in [0,1] (default 0.01) |
| `markovStrandInsensitive` | boolean | No | Table holds strand-insensitive pair frequencies (RSAT 2str; default false) |
| `initialGranularity` | number | No | TFM-Pvalue initial granularity (default 0.1) |
| `maxGranularity` | number | No | TFM-Pvalue finest granularity (default: limited by g·Σmax\|W\| ≤ 2^40) |
| `decreaseFactor` | number | No | Granularity decrease factor (> 1, default 10) |
| `maxStates` | integer | No | DP states per column before a scale is abandoned (default 2097152) |
| `maxSuffixSet` | integer | No | Suffix-score set size for pruning (default 1048576) |
| `exhaustive` | boolean | No | Unbounded enumeration fallback — always exact (default false) |

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
| 1003 | Give either background or markovFrequencies, not both |
| 1003 | markovFrequencies must not be empty |
| 1003 | Invalid search options |
| 1003 | PWM p-values need an explicit word distribution (Markov order ≤ 10) |

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

### Example 3: Dinucleotide (order-1 Markov) background

**Tool Call:**
```json
{
  "tool": "pwm_score_pvalue",
  "arguments": {
    "pwm": { "matrix": [[1.5, -1, 0.75, -0.5], [-0.5, 2, -1.5, 0.5], [0.25, 0.5, 1, -0.75], [-1, -0.25, 0.125, 1.25]], "length": 4 },
    "score": 3.0,
    "markovFrequencies": { "AA": 2, "AC": 1, "AG": 2, "AT": 4, "CA": 2, "CC": 1, "CG": 3, "CT": 2, "GA": 2, "GC": 4, "GG": 2, "GT": 1, "TA": 3, "TC": 2, "TG": 2, "TT": 2 },
    "markovPseudoFrequency": 0
  }
}
```

**Response:**
```json
{
  "score": 3,
  "scoreAboveMaximum": false,
  "pValue": 0.12487874779541444,
  "pValueLowerBound": 0.12487874779541444,
  "pValueUpperBound": 0.12487874779541444,
  "isExact": true,
  "granularity": 10
}
```

MACRO-APE 3.0.6 (`ru.autosome.ape.di.FindPvalue … 3 --from-mono -d 16 -b <the 16 frequencies>`) gives 0.12487874779541447.

### Example 4: Exhaustive mode with a one-state budget

**Tool Call:**
```json
{
  "tool": "pwm_score_pvalue",
  "arguments": {
    "pwm": { "matrix": [[2.0, -2.0], [-2.0, 2.0], [-2.0, -2.0], [-2.0, -2.0]], "length": 2 },
    "score": 0.0, "maxStates": 1, "maxSuffixSet": 1, "exhaustive": true
  }
}
```

**Response:** `{"score": 0, "scoreAboveMaximum": false, "pValue": 0.4375, "pValueLowerBound": 0.4375, "pValueUpperBound": 0.4375, "isExact": true, "granularity": 10}` (without `exhaustive`: bounds 0 … 1, `isExact` false).

## Worked Example

Words AC score 4 (1/16), A[AGT] and [CGT]C score 0 (6/16), all others −4: P(S ≥ 0) = 7/16 = 0.4375. For p = 0.1 the smallest qualifying word score is 4 (P = 1/16); 0 would give 7/16 > 0.1.

## See Also

- [pwm_score_thresholds](pwm_score_thresholds.md) — Biopython grid-approximation thresholds
- [scan_with_pwm](scan_with_pwm.md) — Scan with the chosen threshold
- [create_pwm](create_pwm.md) — Build a PWM
- [alphabet_pwm_score_pvalue](alphabet_pwm_score_pvalue.md) — Same engine for protein / any-alphabet PWMs

## References

- Touzet H., Varré J.-S. (2007). Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2:15. Reference C++ source: CRAN `TFMPvalue` `src/Matrix.cpp`, `src/TFMpvalue.cpp`.
- Vorontsov I.E., Kulakovskiy I.V., Makeev V.J. (2013). Jaccard index based similarity measure to compare transcription factor binding site models. Algorithms Mol Biol 8:23 (MACRO-APE; dinucleotide background `DiPWMScoresGenerator`).
- Boeva V. et al. (2007). Exact p-value calculation for heterotypic clusters of regulatory motifs (AhoPro). Algorithms Mol Biol 2:13 — Markov background p-values.
- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L45](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs#L45)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
