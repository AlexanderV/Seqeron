# integrated_haplotype_score

Compute the integrated haplotype score (iHS) from precomputed EHH curves.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `integrated_haplotype_score` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions)` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Integrates the ancestral (0) and derived (1) EHH decay curves against genomic positions using the trapezoidal rule to obtain iHH0 and iHH1, then returns `iHS = ln(iHH1 / iHH0)`. This overload takes precomputed per-allele EHH values (one per position); it returns 0 when the three arrays are misaligned, hold fewer than two points, or an integrated area is non-positive.

A positive iHS indicates the derived allele has the longer, more conserved haplotype; a negative iHS indicates the ancestral allele does.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L855](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L855)
- Algorithm doc: [Integrated_Haplotype_Score.md](../../../algorithms/Population_Genetics/Integrated_Haplotype_Score.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `ehh0` | number[] | Yes | EHH values for the ancestral (0) allele, aligned with positions |
| `ehh1` | number[] | Yes | EHH values for the derived (1) allele, aligned with positions |
| `positions` | integer[] | Yes | Genomic positions for the EHH entries |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ihs` | number | ln(iHH1 / iHH0); 0 for misaligned/short input or non-positive area |

## Errors

None (misaligned or too-short input returns 0).

## Examples

### Example 1: Extended derived allele → positive iHS

**User Prompt:**
> iHS when the derived allele's EHH stays at 1.0 while the ancestral drops to 0.5 over 1000 bp.

**Expected Tool Call:**
```json
{
  "tool": "integrated_haplotype_score",
  "arguments": { "ehh0": [1.0, 0.5], "ehh1": [1.0, 1.0], "positions": [0, 1000] }
}
```

**Response:**
```json
{ "ihs": 0.287682 }
```

iHH0 = ½(1+0.5)·1000 = 750, iHH1 = ½(1+1)·1000 = 1000; iHS = ln(1000/750) = ln(4/3) ≈ 0.2877.

### Example 2: Balanced EHH decay → iHS = 0

**User Prompt:**
> iHS when both alleles decay identically.

**Expected Tool Call:**
```json
{
  "tool": "integrated_haplotype_score",
  "arguments": {
    "ehh0": [1.0, 0.8, 0.5, 0.2, 0.1],
    "ehh1": [1.0, 0.8, 0.5, 0.2, 0.1],
    "positions": [0, 1000, 2000, 3000, 4000]
  }
}
```

**Response:**
```json
{ "ihs": 0.0 }
```

Identical curves give iHH0 = iHH1, so ln(iHH1/iHH0) = 0.

## Performance

- **Time Complexity:** O(m) in the number of positions
- **Space Complexity:** O(1)

## References

- Voight, B. F., Kudaravalli, S., Wen, X., Pritchard, J. K. (2006). A map of recent positive selection in the human genome. *PLoS Biology* 4(3):e72.

## See Also

- [scan_selection_signals](scan_selection_signals.md)
