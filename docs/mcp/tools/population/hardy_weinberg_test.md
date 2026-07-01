# hardy_weinberg_test

Test Hardy-Weinberg equilibrium for a biallelic variant.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `hardy_weinberg_test` |
| **Method ID** | `PopulationGeneticsAnalyzer.TestHardyWeinberg` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

From observed genotype counts AA/Aa/aa, estimates the allele frequency `p = (2·AA + Aa)/(2n)`, computes the HWE-expected counts `p²n`, `2pqn`, `q²n`, and performs a chi-square goodness-of-fit test with 1 degree of freedom. The variant is reported `inEquilibrium` when the p-value is ≥ the significance level (default 0.05). A zero-sample input returns χ² = 0, p-value = 1, in equilibrium.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L424](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L424)
- Algorithm doc: [Hardy_Weinberg_Test.md](../../../algorithms/Population_Genetics/Hardy_Weinberg_Test.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variantId` | string | Yes | Variant identifier (echoed in result) |
| `observedAA` | integer | Yes | Observed homozygous-AA count |
| `observedAa` | integer | Yes | Observed heterozygous-Aa count |
| `observedaa` | integer | Yes | Observed homozygous-aa count |
| `significanceLevel` | number | No | Significance level α (default 0.05) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variantId` | string | Echoed variant id |
| `observedAA` / `observedAa` / `observedaa` | integer | Echoed observed counts |
| `expectedAA` / `expectedAa` / `expectedaa` | number | HWE-expected counts (p²n, 2pqn, q²n) |
| `chiSquare` | number | Chi-square statistic |
| `pValue` | number | P-value (1 df) |
| `inEquilibrium` | boolean | True iff pValue ≥ significanceLevel |

## Errors

None (zero samples returns χ² = 0, p-value = 1).

## Examples

### Example 1: Ford's scarlet tiger moth data (in equilibrium)

**User Prompt:**
> Test HWE for Ford's moth counts: 1469 AA, 138 Aa, 5 aa.

**Expected Tool Call:**
```json
{
  "tool": "hardy_weinberg_test",
  "arguments": { "variantId": "FORD_MOTH", "observedAA": 1469, "observedAa": 138, "observedaa": 5 }
}
```

**Response (key fields):**
```json
{
  "expectedAA": 1467.40,
  "expectedAa": 141.21,
  "expectedaa": 3.40,
  "chiSquare": 0.8309,
  "inEquilibrium": true
}
```

χ² ≈ 0.83 < 3.841 (α = 0.05, 1 df) → population is in Hardy-Weinberg equilibrium.

### Example 2: Excess heterozygotes (deviates from HWE)

**User Prompt:**
> Test HWE for 10 AA, 80 Aa, 10 aa.

**Expected Tool Call:**
```json
{
  "tool": "hardy_weinberg_test",
  "arguments": { "variantId": "EXCESS_HET", "observedAA": 10, "observedAa": 80, "observedaa": 10 }
}
```

**Response (key fields):**
```json
{
  "expectedAA": 25.0,
  "expectedAa": 50.0,
  "expectedaa": 25.0,
  "chiSquare": 36.0,
  "inEquilibrium": false
}
```

χ² = (10−25)²/25 + (80−50)²/50 + (10−25)²/25 = 9 + 18 + 9 = 36 → rejects HWE.

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## References

- Hardy, G. H. (1908); Weinberg, W. (1908).
- Ford, E. B. (1971). *Ecological Genetics* (scarlet tiger moth dataset).
- Wikipedia contributors. [Hardy–Weinberg principle](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle).

## See Also

- [allele_frequencies](allele_frequencies.md)
