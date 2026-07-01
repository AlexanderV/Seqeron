# epigenetic_age

Estimate epigenetic age (Horvath-clock-style) from methylation at clock CpGs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `epigenetic_age` |
| **Method ID** | `EpigeneticsAnalyzer.CalculateEpigeneticAge` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes a Horvath (2013) epigenetic-clock age estimate. The linear predictor is
`intercept + Σ coef_i · β_i` over the CpGs shared by the supplied methylation β-values and coefficient
table (CpGs absent from either side contribute nothing). The predictor is mapped to years with the Horvath
inverse calibration `F⁻¹` (`adult.age = 20`):

```
antiTransform(x) = (1 + 20)·exp(x) − 1        for x < 0
antiTransform(x) = (1 + 20)·x + 20            for x ≥ 0
```

The caller provides the published clock coefficient table, keeping the wrapper independent of any specific
built-in clock.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L1171](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L1171)
- Inverse transform: [EpigeneticsAnalyzer.cs#L1206](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L1206)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `methylationAtClockCpGs` | object | Yes | CpG ID → methylation β-value (0..1) |
| `coefficients` | object | Yes | CpG ID → clock coefficient |
| `intercept` | number | No | Model intercept (default 0.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `age` | number | Estimated DNAm age in years |

## Errors

| Code | Message |
|------|---------|
| 1001 | Methylation map cannot be null or empty |
| 1001 | Clock coefficient table cannot be null or empty |

## Examples

### Example 1: Positive linear predictor

With intercept 0.5 and predictor 0.5 + 1.0·0.5 + 2.0·0.25 = 1.5:

**Response:**
```json
{ "age": 51.5 }
```

`(1 + 20)·1.5 + 20 = 51.5`.

### Example 2: Negative predictor

Predictor −1.0 uses the exponential branch:

**Response:**
```json
{ "age": 6.7254682646002895 }
```

`(1 + 20)·exp(−1) − 1 ≈ 6.7255`.

## Performance

- **Time Complexity:** O(k) for k supplied CpGs
- **Space Complexity:** O(1)

## See Also

- [methylation_profile](methylation_profile.md) — aggregate methylation levels
- [methylation_from_bisulfite](methylation_from_bisulfite.md) — per-CpG methylation from reads
