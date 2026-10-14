# thermodynamics

Calculate thermodynamic properties of a DNA duplex.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `thermodynamics` |
| **Method ID** | `SequenceStatistics.CalculateThermodynamics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates thermodynamic properties of a DNA duplex using the nearest-neighbor method. Returns enthalpy (ΔH), entropy (ΔS), Gibbs free energy (ΔG), and melting temperature (Tm). Reproduces Biopython 1.88 `Tm_NN(nn_table=DNA_NN3)` (dnac1 = dnac2 = C_T/2; self-complementary: dnac1 = C_T). Default: SantaLucia (1998) method-5 Na+ entropy correction, rounded output. Optional: self-complementary duplex, any Biopython `salt_correction` method 0–7 with K+/Tris/Mg2+/dNTPs (all in mol/L, like `naConcentration`), unrounded output. For methods other than 5, `deltaS` is the 1 M NaCl value (the correction acts on Tm or 1/Tm).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L1329](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L1329)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence (min length: 1) |
| `naConcentration` | number | No | Na+ concentration in M (default: 0.05 = 50mM) |
| `primerConcentration` | number | No | Total strand concentration C_T in M (default: 0.00000025 = 250 nM — a repository convention; Biopython `Tm_NN` defaults to 50 nM total) |
| `selfComplementary` | boolean | No | Self-complementary (homo)duplex: ΔS° −1.4 symmetry correction, Tm uses C_T/1 (Biopython `selfcomp=True`; default false) |
| `saltMethod` | integer | No | Biopython salt-correction method 0–7 (0 none, 1–4 on Tm, 5 SantaLucia 1998 entropy = default, 6 Owczarzy 2004, 7 Owczarzy 2008 Mg²⁺/dNTPs) |
| `kConcentration` | number | No | K+ concentration in M (default 0) |
| `trisConcentration` | number | No | Tris concentration in M (default 0) |
| `mgConcentration` | number | No | Mg2+ concentration in M (default 0) |
| `dntpsConcentration` | number | No | dNTPs concentration in M (default 0) |
| `roundResults` | boolean | No | Round ΔH/ΔS/ΔG to 2 dp and Tm to 1 dp (default true); false = unrounded |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `deltaH` | number | Enthalpy change (kcal/mol) |
| `deltaS` | number | Entropy change (cal/mol·K) |
| `deltaG` | number | Gibbs free energy change at 37°C (kcal/mol) |
| `meltingTemperature` | number | Melting temperature (°C) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Neighbour pair without nearest-neighbour parameters (inosine, e.g. `II`; Biopython `Tm_NN` raises ValueError) |
| 1001 | Concentration not positive (C_T, default-method Na+) / negative or non-finite, `saltMethod` outside 0–7, or zero total ions for methods 1–6 |

## Examples

### Example 1: DNA duplex thermodynamics

**User Prompt:**
> What are the thermodynamic properties of "ATGCGATCGATCG"?

**Expected Tool Call:**
```json
{
  "tool": "thermodynamics",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "deltaH": -102.1,
  "deltaS": -290.63,
  "deltaG": -11.96,
  "meltingTemperature": 42.4
}
```

### Example 2: Custom salt concentration

**User Prompt:**
> Calculate Tm for "GCGCGCGC" at 100mM Na+

**Expected Tool Call:**
```json
{
  "tool": "thermodynamics",
  "arguments": {
    "sequence": "GCGCGCGC",
    "naConcentration": 0.1
  }
}
```

### Example 3: Mg2+/dNTPs (Owczarzy 2008)

**Expected Tool Call:**
```json
{
  "tool": "thermodynamics",
  "arguments": {
    "sequence": "CGTTCCAAAGATGTGGGCATGAGCTTAC",
    "saltMethod": 7,
    "kConcentration": 0.02,
    "trisConcentration": 0.01,
    "mgConcentration": 0.0015,
    "dntpsConcentration": 0.0008
  }
}
```

**Response:** `meltingTemperature` 67.7 (Biopython 1.88 `Tm_NN(..., Na=50, K=20, Tris=10, Mg=1.5, dNTPs=0.8, dnac1=dnac2=125, saltcorr=7)` = 67.73409261682195; `roundResults: false` returns that value).

### Example 4: Self-complementary duplex

`{"sequence": "CGCGAATTCGCG", "selfComplementary": true}` → `meltingTemperature` 47.0 (Biopython `selfcomp=True`, dnac1 = 250 nM: 46.96872419752236).

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [melting_temperature](melting_temperature.md) - Simple Tm calculation
- [nucleotide_composition](nucleotide_composition.md) - GC content analysis
