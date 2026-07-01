# analyze_oligo

Analyze the basic physical properties of a short oligonucleotide.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `analyze_oligo` |
| **Method ID** | `ProbeDesigner.AnalyzeOligo` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the melting temperature (Tm), GC content, molecular weight, and 260 nm molar
extinction coefficient of an oligonucleotide (primer, probe, or short synthetic oligo).
Call this when the user needs the basic physical characterization of an oligo.

- **Tm** uses the **Wallace rule** (`2·(A+T) + 4·(G+C)`) for sequences shorter than 14 bases,
  and a **salt-adjusted formula** (`81.5 + 16.6·log10(0.05) + 41·%GC − 600/length`, 50 mM Na⁺)
  for sequences of 14 bases or longer.
- **GC content** is returned as a **fraction in [0,1]** (not a percentage).
- **Molecular weight** (Da) sums per-base average weights (A=331.2, C=307.2, G=347.2,
  T=322.2, U=308.2) and subtracts one water (18.0 Da) per phosphodiester bond.
- **Extinction coefficient** (M⁻¹·cm⁻¹ at 260 nm) sums per-base contributions
  (A=15400, C=7400, G=11500, T=8700, U=9900; unknown bases fall back to 10000).

Input is treated case-insensitively.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1299](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1299)
- Tm constants: [ThermoConstants.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Oligonucleotide sequence (A/C/G/T/U, case-insensitive; min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tm` | number | Melting temperature in °C |
| `gcContent` | number | GC content as a fraction (0–1) |
| `molecularWeight` | number | Molecular weight in Da |
| `extinctionCoefficient` | number | Molar extinction coefficient at 260 nm (M⁻¹·cm⁻¹) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: 4-mer (Wallace rule Tm)

**User Prompt:**
> Analyze the oligo "ATGC".

**Expected Tool Call:**
```json
{
  "tool": "analyze_oligo",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "tm": 12,
  "gcContent": 0.5,
  "molecularWeight": 1253.8,
  "extinctionCoefficient": 43000
}
```

### Example 2: 20-mer (salt-adjusted Tm)

**User Prompt:**
> What are the properties of "ACGTACGTACGTACGTACGT"?

**Expected Tool Call:**
```json
{
  "tool": "analyze_oligo",
  "arguments": {
    "sequence": "ACGTACGTACGTACGTACGT"
  }
}
```

**Response:**
```json
{
  "tm": 50.402902071977906,
  "gcContent": 0.5,
  "molecularWeight": 6197,
  "extinctionCoefficient": 215000
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [oligo_extinction_coefficient](oligo_extinction_coefficient.md) - Extinction coefficient only
- [oligo_concentration_from_absorbance](oligo_concentration_from_absorbance.md) - Beer–Lambert concentration
- [primer_melting_temperature](primer_melting_temperature.md) - Standalone Tm calculation
