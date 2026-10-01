# create_alphabet_pwm

Log-odds PWM over any alphabet (protein, RNA, gapped DNA).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `create_alphabet_pwm` |
| **Method ID** | `MotifFinder.CreateAlphabetPwm` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds a log-odds position weight matrix over an arbitrary alphabet from aligned instances — Biopython `motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)`: W[a,j] = log2(((c[a,j] + p[a]) / (Σ c[·,j] + Σ p)) / q[a]). Returns the K×L matrix (rows in alphabet order), the consensus / anticonsensus (first symbol in alphabet order on ties, as Biopython), the maximal / minimal score (Biopython `max` / `min`) and the mean / standard deviation of the score of a random background window (Biopython `mean` / `std`, −∞ cells skipped). Matching is case-insensitive. Symbols outside the alphabet are rejected, or skipped when `ignoreUnknownSymbols` (Biopython keeps only alphabet letters). JSON has no −∞: matrix cells and `minScore` are null where the value is −∞ (unseen symbol with zero pseudocount).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L42](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L42)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned instances of equal length |
| `alphabet` | string | Yes | Alphabet in row order, distinct ignoring case (e.g. ACDEFGHIKLMNPQRSTVWY) |
| `pseudocount` | number | No | Pseudocount per cell (default 0 = Biopython None); ignored when pseudocounts is given |
| `pseudocounts` | array<number> | No | Per-symbol pseudocounts in alphabet order |
| `background` | array<number> | No | Background in alphabet order (> 0, normalised; uniform when omitted); also used for mean/std |
| `ignoreUnknownSymbols` | boolean | No | Skip symbols outside the alphabet when counting (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alphabet` | string | Alphabet (row order) |
| `matrix` | array<array<number|null>> | K×L log-odds; null = −∞ |
| `length` | integer | Motif length L |
| `consensus` | string | Consensus |
| `anticonsensus` | string | Anticonsensus |
| `maxScore` | number|null | Σ column maxima |
| `minScore` | number|null | Σ column minima; null = −∞ |
| `mean` | number | Expected background score |
| `std` | number | Standard deviation of the background score |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required. |
| 1001 | Alphabet is required (e.g. ACDEFGHIKLMNPQRSTVWY). |
| 1001 | Invalid character ... not in the alphabet |
| 1001 | Alphabet symbol ... occurs more than once (symbols are compared ignoring case). |
| 1002 | Pseudocount must be finite and non-negative. |

## Examples

### Example 1: Three short protein instances

**User Prompt:**
> Build a PWM from MKV, MRV, LKI over the alphabet IKLMRV with pseudocount 0.5.

**Tool Call:**
```json
{
  "tool": "create_alphabet_pwm",
  "arguments": {
    "sequences": [
      "MKV",
      "MRV",
      "LKI"
    ],
    "alphabet": "IKLMRV",
    "pseudocount": 0.5
  }
}
```

**Response:**
```json
{
  "alphabet": "IKLMRV",
  "matrix": [
    [
      -1,
      -1,
      0.5849625007211562
    ],
    [
      -1,
      1.3219280948873626,
      -1
    ],
    [
      0.5849625007211562,
      -1,
      -1
    ],
    [
      1.3219280948873626,
      -1,
      -1
    ],
    [
      -1,
      0.5849625007211562,
      -1
    ],
    [
      -1,
      -1,
      1.3219280948873626
    ]
  ],
  "length": 3,
  "consensus": "MKV",
  "anticonsensus": "IIK",
  "maxScore": 3.965784284662088,
  "minScore": -3,
  "mean": 1.0911319941500706,
  "std": 1.7447483665609653
}
```

## Worked Example

Every field equals Biopython 1.88 (`motifs.create(['MKV','MRV','LKI'], alphabet='IKLMRV').counts.normalize(0.5).log_odds()`: consensus MKV, anticonsensus IIK, max 3.965784284662088, min −3, mean 1.0911319941500706, std 1.7447483665609653). 400 random protein / sub-alphabet / DNA cases (scalar and per-symbol pseudocounts, backgrounds, ignored X/− symbols) agree with Biopython to ≤ 1.1e-14 relative.

## See Also

- [scan_with_alphabet_pwm](scan_with_alphabet_pwm.md) — Score / scan a sequence with such a PWM
- [create_pwm](create_pwm.md) — DNA (A,C,G,T) PWM

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L42](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs#L42)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
