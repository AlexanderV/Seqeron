# scan_selection_signals

Scan genomic regions for selection signals using Tajima's D, Fst, and iHS thresholds.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `scan_selection_signals` |
| **Method ID** | `PopulationGeneticsAnalyzer.ScanForSelection(regions, thresholds)` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each region with precomputed Tajima's D, Fst, and iHS scores, emits a signal for every test that crosses its threshold:

- `tajimaD < tajimaDThreshold` → a `TajimasD` signal (excess rare variants),
- `fst > fstThreshold` → an `Fst` signal (high differentiation / local adaptation),
- `|iHS| > ihsThreshold` → an `iHS` signal (haplotype-based selection; direction from the sign of iHS).

A region may yield zero, one, or several signals. Each carries an approximate p-value and an interpretation string.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L882](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L882)
- Algorithm doc: [Integrated_Haplotype_Score.md](../../../algorithms/Population_Genetics/Integrated_Haplotype_Score.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `regions` | object[] | Yes | Regions with `region`, `start`, `end`, `tajimaD`, `fst`, `iHS` |
| `tajimaDThreshold` | number | No | Tajima's D threshold (default -2.0) |
| `fstThreshold` | number | No | Fst threshold (default 0.25) |
| `ihsThreshold` | number | No | |iHS| threshold (default 2.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | object[] | Emitted signals |
| `items[].region` | string | Region id |
| `items[].start` / `items[].end` | integer | Region bounds |
| `items[].score` | number | The score that crossed the threshold |
| `items[].testType` | string | `TajimasD`, `Fst`, or `iHS` |
| `items[].pValue` | number | Approximate p-value |
| `items[].interpretation` | string | Human-readable interpretation |

## Errors

None (a region crossing no threshold simply emits nothing).

## Examples

### Example 1: Strongly negative Tajima's D emits a signal

**User Prompt:**
> Scan a region with Tajima's D = -2.5.

**Expected Tool Call:**
```json
{
  "tool": "scan_selection_signals",
  "arguments": {
    "regions": [{ "region": "Region1", "start": 0, "end": 10000, "tajimaD": -2.5, "fst": 0.1, "iHS": 0.5 }],
    "tajimaDThreshold": -2.0
  }
}
```

**Response:** one `TajimasD` signal (score -2.5) with interpretation "Possible positive/purifying selection (excess rare variants)".

### Example 2: Neutral region emits no signals

**User Prompt:**
> Scan a region with Tajima's D 0, Fst 0.1, iHS 0.5.

**Expected Tool Call:**
```json
{
  "tool": "scan_selection_signals",
  "arguments": {
    "regions": [{ "region": "Region1", "start": 0, "end": 10000, "tajimaD": 0.0, "fst": 0.1, "iHS": 0.5 }]
  }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(regions)
- **Space Complexity:** O(signals)

## References

- Voight, B. F. et al. (2006). *PLoS Biology* 4(3):e72.
- Tajima, F. (1989). *Genetics* 123:585–595.

## See Also

- [integrated_haplotype_score](integrated_haplotype_score.md), [tajimas_d](tajimas_d.md), [fst](fst.md)
