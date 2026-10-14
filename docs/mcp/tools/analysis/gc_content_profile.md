# gc_content_profile

GC content in sliding windows along a nucleotide sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `gc_content_profile` |
| **Method ID** | `SequenceStatistics.CalculateGcContentProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **GC-content profile**: the GC content of each sliding window of `windowSize`
bases, advanced by `stepSize` — a percentage (0–100) by default or a fraction (0–1, Biopython
`gc_fraction` scale) with `fraction = true`.

- **Default counting** (no `ambiguityMode`): `(G + C) / (A + C + G + T + U)` over the window; every other
  symbol (N, IUPAC codes, gaps) is excluded from numerator and denominator. Equals Biopython
  `gc_fraction(window, "remove")` for windows without S/W.
- **`ambiguityMode`** (Biopython `gc_fraction(…, ambiguous=…)`, case-insensitive): `remove` (S counts as GC,
  S/W in the denominator, other codes dropped), `ignore` (denominator = window length), `weighted`
  (each code adds its mean GC content, N = 0.5, to the numerator; denominator = window length).
- **`includePartialWindow`**: after the complete windows, every further start `k·stepSize < length` is
  scored on the truncated trailing slice — Biopython `GC_skew`-style `range(0, len(seq), step)` slicing
  (also when `windowSize` exceeds the length). Off (default): only complete windows; a window longer
  than the sequence gives an empty profile.

The windows are produced by the canonical `GcSkewCalculator.CalculateWindowedGcContent` driver.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L2416](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2416) (default), [#L2151](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2443) (ambiguity mode), [#L2191](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2483) / [#L2220](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2512) (partial windows)
- Algorithm: `docs/algorithms/Statistics/GC_Content_Profile.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1) |
| `windowSize` | integer | No | Window size (default 100, ≥ 1) |
| `stepSize` | integer | No | Step size (default 1, ≥ 1) |
| `fraction` | boolean | No | `true`: values in [0, 1]; `false` (default): percentages |
| `ambiguityMode` | string | No | `remove` / `ignore` / `weighted` (Biopython `gc_fraction`); omit for the default counting |
| `includePartialWindow` | boolean | No | `true`: also score trailing windows shorter than `windowSize` (default `false`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | array of number | GC content per window (percentage 0–100, or fraction 0–1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |
| 1002 | Step size must be at least 1 |
| 1003 | Unknown ambiguity mode '…'. Expected one of: remove, ignore, weighted. |

## Examples

### Example 1: GC-rich windows (window 2, step 1)

**User Prompt:**
> GC profile of "GGCC" with window 2.

**Expected Tool Call:**
```json
{
  "tool": "gc_content_profile",
  "arguments": { "sequence": "GGCC", "windowSize": 2, "stepSize": 1 }
}
```

**Response:**
```json
{ "values": [100.0, 100.0, 100.0] }
```
Windows GG, GC, CC are all 100% GC.

### Example 2: AT-only windows, non-overlapping (window 2, step 2)

**User Prompt:**
> GC profile of "ATAT" with window 2, step 2.

**Expected Tool Call:**
```json
{
  "tool": "gc_content_profile",
  "arguments": { "sequence": "ATAT", "windowSize": 2, "stepSize": 2 }
}
```

**Response:**
```json
{ "values": [0.0, 0.0] }
```
Both non-overlapping windows (AT, AT) contain no G or C.

### Example 3: Biopython ambiguity modes as fractions

**Input:** `{ "sequence": "ACGTSSWWNNRYacgusw", "windowSize": 6, "stepSize": 3, "fraction": true, "ambiguityMode": "weighted" }`

**Response:** `{ "values": [0.6666666666666666, 0.4166666666666667, 0.3333333333333333, 0.5833333333333334, 0.5] }`
(Biopython 1.88 `gc_fraction(s[i:i+6], "weighted")`; `remove` → 0.667, 0.4, 0, 0.667, 0.5; `ignore` → 0.667, 0.333, 0, 0.333, 0.5.)

### Example 4: Trailing partial windows

**Input:** `{ "sequence": "GGGCCCAAATTTGCGCA", "windowSize": 6, "stepSize": 4, "fraction": true, "includePartialWindow": true }`

**Response:** `{ "values": [1.0, 0.3333333333333333, 0.3333333333333333, 0.8, 0.0] }`
(Starts 0, 4, 8 are complete; 12 (`GCGCA`) and 16 (`A`) are partial = `[gc_fraction(s[i:i+6]) for i in range(0, 17, 4)]`.)

## Performance

- **Time Complexity:** O(n · windowSize / stepSize).
- **Space Complexity:** O(number of windows).

## See Also

- [gc_content](../sequence/gc_content.md) — whole-sequence GC content
- [analyze_gc_content](analyze_gc_content.md) — comprehensive GC report
- [entropy_profile](entropy_profile.md) — windowed Shannon entropy
