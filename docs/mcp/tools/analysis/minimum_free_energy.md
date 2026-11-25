# minimum_free_energy

Zuker minimum free energy of an RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `minimum_free_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateMinimumFreeEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **minimum free energy** (MFE) of an RNA sequence with a Zuker-style
dynamic program using Turner 2004 nearest-neighbour parameters (O(n³)). The result is
in kcal/mol; a sequence that cannot form any structure returns 0. Loops shorter than
`minLoopSize` (default 3) are disallowed; smaller values are clamped to 3.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L1020](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L1020)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (min length 1) |
| `minLoopSize` | integer | No | Minimum hairpin loop size (default 3; values < 3 clamped) |
| `temperature` | number | No | Absolute temperature in Kelvin (default 310.15 = 37 °C); every Turner 2004 parameter is rescaled with its enthalpy, dG(T) = dH - (dH - dG37)*T/310.15 (ViennaRNA md.temperature); must be positive |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `mfe` | number | Minimum free energy in kcal/mol |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Simple GC hairpin

**User Prompt:**
> What is the MFE of "GGGAAACCC"?

**Expected Tool Call:**
```json
{
  "tool": "minimum_free_energy",
  "arguments": { "rnaSequence": "GGGAAACCC" }
}
```

**Response:**
```json
{ "mfe": -1.12 }
```
Three GC pairs (stacking −6.52) plus a 3-nt hairpin initiation (+5.4) ⇒ −1.12 kcal/mol.

### Example 2: No structure (poly-A)

**User Prompt:**
> MFE of "AAAAAAAA".

**Expected Tool Call:**
```json
{
  "tool": "minimum_free_energy",
  "arguments": { "rnaSequence": "AAAAAAAA" }
}
```

**Response:**
```json
{ "mfe": 0.0 }
```

## Performance

- **Time Complexity:** O(n³) (multiloop split) + O(n²·30²) (interior loops ≤ 30 nt); ~0.8 s at n = 500, ~3.5 s at n = 1000 (Release).
- **Space Complexity:** O(n²).

## See Also

- [rna_mfe_structure](rna_mfe_structure.md) — the MFE structure (dot-bracket + pairs) for this value
- [predict_rna_structure](predict_rna_structure.md) — stem-loop assembly (not the MFE fold)
- [find_stem_loops](find_stem_loops.md) — hairpin enumeration with energies
