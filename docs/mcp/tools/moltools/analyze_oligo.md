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

- **Tm** is Primer3's `seqtm` (identical to primer3-py `calc_tm(..., max_nn_length=36)`, Primer3 MAX_NN_TM_LENGTH) at the Primer3
  hybridization-probe defaults (50 nM oligo, 50 mM monovalent, no Mg²⁺, no dNTP): SantaLucia 1998
  nearest-neighbour Tm with the SantaLucia salt correction for ≤ 36 nt, Primer3 `long_seq_tm`
  (`81.5 + 16.6·log10([Na⁺]) + 0.41·%GC − 600/N`) above. It is `null` when not computable
  (fewer than 2 bases, or any non-ACGT base — e.g. an RNA oligo).
- **GC content** is returned as a **fraction in [0,1]** (not a percentage).
- **Molecular weight** (Da) is the single-stranded Biopython `molecular_weight(seq, "DNA"|"RNA")`
  (canonical `SequenceStatistics.CalculateNucleotideMolecularWeight`; nucleoside-monophosphate masses
  minus one water 18.01528 Da per phosphodiester bond; RNA when the oligo contains U and no T;
  unknown symbols skipped).
- **Extinction coefficient** (M⁻¹·cm⁻¹ at 260 nm) sums per-base contributions
  (A=15400, C=7400, G=11500, T=8700, U=9900; unknown bases fall back to 10000).

Input is treated case-insensitively.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L3506](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L3506) (`AnalyzeOligo`)
- Tm constants: [ThermoConstants.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Oligonucleotide sequence (A/C/G/T/U, case-insensitive; min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tm` | number \| null | Melting temperature in °C (null when not computable) |
| `gcContent` | number | GC content as a fraction (0–1) |
| `molecularWeight` | number | Molecular weight in Da |
| `extinctionCoefficient` | number | Molar extinction coefficient at 260 nm (M⁻¹·cm⁻¹) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: 4-mer

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
  "tm": -54.53620887888866,
  "gcContent": 0.5,
  "molecularWeight": 1253.8027,
  "extinctionCoefficient": 43000
}
```

### Example 2: 20-mer (nearest-neighbour Tm)

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
  "tm": 53.99351583691026,
  "gcContent": 0.5,
  "molecularWeight": 6196.9523,
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
