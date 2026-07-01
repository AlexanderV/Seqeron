# analyze_telomeres

Analyze 5'/3' telomere repeat tracts on a chromosome sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `analyze_telomeres` |
| **Method ID** | `ChromosomeAnalyzer.AnalyzeTelomeres` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans the two chromosome ends for telomeric repeats. The 3' end is matched against the repeat unit
(`TTAGGG` by default); the 5' end is matched against its reverse complement (`CCCTAA`). Repeat
windows are walked inward from each end while per-window similarity to the unit stays ≥ 0.7, summing
the tract length and repeat purity. A telomere is "present" when its length ≥ `minTelomereLength`,
and "critically short" when a present telomere is shorter than `criticalLength`.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L393](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L393)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosomeName` | string | Yes | Chromosome name (echoed in result). |
| `sequence` | string | Yes | Chromosome nucleotide sequence. |
| `telomereRepeat` | string | No | Repeat unit (default `TTAGGG`). |
| `searchLength` | integer | No | End-region search length in bp (default 10000, > 0). |
| `minTelomereLength` | integer | No | Minimum length to call a telomere present (default 500). |
| `criticalLength` | integer | No | Critical shortening threshold (default 3000). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | string | Echoed name. |
| `has5PrimeTelomere` / `has3PrimeTelomere` | boolean | Telomere present at each end. |
| `telomereLength5Prime` / `telomereLength3Prime` | integer | Tract length in bp. |
| `repeatPurity5Prime` / `repeatPurity3Prime` | number | Fraction of matching bases in the tract. |
| `isCriticallyShort` | boolean | True when a present telomere is below `criticalLength`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Chromosome name cannot be null or empty |
| 1003 | Telomere repeat cannot be null or empty |
| 1004 | Search length must be positive |

## Example

100 tandem copies of `TTAGGG` (600 bp): the 3' end is a pure telomere; the 5' end does not start
with `CCCTAA` so no 5' telomere is called.

```json
{
  "chromosome": "chrT",
  "has5PrimeTelomere": false,
  "telomereLength5Prime": 0,
  "has3PrimeTelomere": true,
  "telomereLength3Prime": 600,
  "repeatPurity3Prime": 1.0,
  "isCriticallyShort": true
}
```

## References

- [ChromosomeAnalyzer.AnalyzeTelomeres](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L393)
