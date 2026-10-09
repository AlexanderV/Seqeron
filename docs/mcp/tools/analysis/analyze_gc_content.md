# analyze_gc_content

Comprehensive GC analysis of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `analyze_gc_content` |
| **Method ID** | `GcSkewCalculator.AnalyzeGcContent` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces a comprehensive GC report for a DNA sequence: the overall GC content
(`(G+C)/(A+T+G+C)·100`), overall GC skew (`(G−C)/(G+C)`), overall AT skew
(`(A−T)/(A+T)`), the **population** variance (`Σ(xᵢ−μ)²/N`) of the per-window GC
content and GC skew, and the full sliding-window GC-skew and GC-content profiles.
GC skew is used to locate replication origins/termini in bacterial genomes; the
windowed profiles reveal compositional heterogeneity along the sequence.

Only `A`, `C`, `G`, `T` are counted (case-insensitive); other symbols are ignored
in both numerator and denominator. When the sequence is shorter than `windowSize`
no full window exists, so both windowed profiles are empty and both window-derived
variances are 0, while the overall scalar metrics are still computed over the whole
sequence.

With `fraction=true` (Core `AnalyzeGcContent(…, fraction: true)`) the overall and windowed GC content
are reported in [0,1] like Biopython 1.88 `Bio.SeqUtils.gc_fraction`, and `gcContentVariance` is the
variance of those fractions (= the percentage variance / 10⁴); the skew fields are unchanged.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L861](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L861)
- Evidence: `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T, min length 1) |
| `windowSize` | integer | No | Sliding-window length for the profiles (default 1000) |
| `stepSize` | integer | No | Step between window starts (default 100) |
| `fraction` | boolean | No | Report overall and windowed GC content (and its variance) as a fraction in [0,1] (Biopython gc_fraction) instead of a percentage. Default false. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `overallGcContent` | number | GC% over the whole sequence, `(G+C)/(A+T+G+C)·100` (a fraction in [0,1] with `fraction=true`) |
| `overallGcSkew` | number | GC skew `(G−C)/(G+C)`, in [−1, 1] |
| `overallAtSkew` | number | AT skew `(A−T)/(A+T)`, in [−1, 1] |
| `gcContentVariance` | number | Population variance of per-window GC% (of per-window fractions with `fraction=true`) |
| `gcSkewVariance` | number | Population variance of per-window GC skew |
| `windowedGcSkew` | array | Per-window GC-skew points (`position`, `gcSkew`, `windowStart`, `windowEnd`) |
| `windowedGcContent` | array | Per-window GC-content points (`position`, `gcContent`, `windowStart`, `windowEnd`) |
| `sequenceLength` | integer | Length of the input sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid DNA sequence |

## Examples

### Example 1: Overall scalar metrics

**User Prompt:**
> Give me a full GC analysis of "GGGCCAT".

**Expected Tool Call:**
```json
{
  "tool": "analyze_gc_content",
  "arguments": { "sequence": "GGGCCAT" }
}
```

**Response (scalar fields; profiles empty because len < 1000):**
```json
{
  "overallGcContent": 71.42857142857143,
  "overallGcSkew": 0.2,
  "overallAtSkew": 0.0,
  "gcContentVariance": 0.0,
  "gcSkewVariance": 0.0,
  "windowedGcSkew": [],
  "windowedGcContent": [],
  "sequenceLength": 7
}
```
`GC% = (3+2)/7·100 = 71.428…`, `GC skew = (3−2)/(3+2) = 0.2`, `AT skew = (1−1)/2 = 0`.

### Example 2: Windowed population variance

**User Prompt:**
> Analyze "GGCC" with window 2, step 2.

**Expected Tool Call:**
```json
{
  "tool": "analyze_gc_content",
  "arguments": { "sequence": "GGCC", "windowSize": 2, "stepSize": 2 }
}
```

**Response (key fields):**
```json
{
  "overallGcContent": 100.0,
  "overallGcSkew": 0.0,
  "gcSkewVariance": 1.0,
  "gcContentVariance": 0.0,
  "sequenceLength": 4
}
```
Windows `GG` (skew +1) and `CC` (skew −1): population variance of {+1, −1} is
`((1−0)²+(−1−0)²)/2 = 1.0` (division by N, not N−1). Both windows are 100% GC, so
`gcContentVariance = 0`.

### Example 3: GC content as a fraction

**Expected Tool Call:**
```json
{
  "tool": "analyze_gc_content",
  "arguments": { "sequence": "GGGCCAT", "windowSize": 4, "stepSize": 3, "fraction": true }
}
```

**Response (key fields):**
```json
{
  "overallGcContent": 0.7142857142857143,
  "overallGcSkew": 0.2,
  "gcContentVariance": 0.0625,
  "gcSkewVariance": 0.5625,
  "sequenceLength": 7
}
```
Windows `GGGC` and `CCAT`: Biopython `gc_fraction` 1.0 and 0.5 (population variance 0.0625),
`GC_skew` 0.5 and −1.0 (population variance 0.5625).

## Performance

- **Time Complexity:** O(n · windowSize / stepSize) for the windowed profiles.
- **Space Complexity:** O(number of windows).

## See Also

- [gc_skew](gc_skew.md) — overall GC skew only
- [windowed_gc_skew](windowed_gc_skew.md) — sliding-window GC skew profile
- [predict_replication_origin](predict_replication_origin.md) — origin/terminus from cumulative skew
