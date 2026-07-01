# linkage_disequilibrium

Compute pairwise linkage disequilibrium (D' and r²) between two variants.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `linkage_disequilibrium` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateLD` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

From per-individual genotype pairs (0/1/2 encoding), computes:

- **r²** — the squared Pearson correlation of the genotype values (diploid correlation equals the haplotype correlation; Hill & Robertson 1968), clamped to [0, 1].
- **D'** — from the genotype covariance `D = Cov(X₁,X₂)/2`, normalized by `D_max` and clamped to [0, 1] (Lewontin 1964).

Identical genotype vectors give r² = 1 and D' = 1 (perfect LD); an independent (balanced) design gives r² = 0 and D' = 0. Monomorphic loci and empty input return 0 without NaN. Variant ids and distance are echoed unchanged.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L729](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L729)
- Algorithm doc: [Linkage_Disequilibrium.md](../../../algorithms/Population_Genetics/Linkage_Disequilibrium.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant1Id` | string | Yes | Variant 1 identifier |
| `variant2Id` | string | Yes | Variant 2 identifier |
| `genotypes` | object[] | Yes | Per-individual `{geno1, geno2}` pairs (0/1/2) |
| `distance` | integer | Yes | Genomic distance between variants |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variant1` / `variant2` | string | Echoed variant ids |
| `dPrime` | number | Normalized |D'| ∈ [0, 1] |
| `rSquared` | number | Squared genotype correlation r² ∈ [0, 1] |
| `distance` | number | Echoed genomic distance |

## Errors

None (monomorphic/empty input returns 0).

## Examples

### Example 1: Perfect LD

**User Prompt:**
> LD between two variants with identical genotype vectors.

**Expected Tool Call:**
```json
{
  "tool": "linkage_disequilibrium",
  "arguments": {
    "variant1Id": "V1", "variant2Id": "V2",
    "genotypes": [
      { "geno1": 0, "geno2": 0 }, { "geno1": 0, "geno2": 0 },
      { "geno1": 1, "geno2": 1 }, { "geno1": 1, "geno2": 1 },
      { "geno1": 2, "geno2": 2 }, { "geno1": 2, "geno2": 2 }
    ],
    "distance": 1000
  }
}
```

**Response:**
```json
{ "variant1": "V1", "variant2": "V2", "dPrime": 1.0, "rSquared": 1.0, "distance": 1000 }
```

### Example 2: No LD (balanced design)

**User Prompt:**
> LD for two variants where every genotype combination appears equally.

**Expected Tool Call:**
```json
{
  "tool": "linkage_disequilibrium",
  "arguments": {
    "variant1Id": "V1", "variant2Id": "V2",
    "genotypes": [
      { "geno1": 0, "geno2": 0 }, { "geno1": 0, "geno2": 1 }, { "geno1": 0, "geno2": 2 },
      { "geno1": 1, "geno2": 0 }, { "geno1": 1, "geno2": 1 }, { "geno1": 1, "geno2": 2 },
      { "geno1": 2, "geno2": 0 }, { "geno1": 2, "geno2": 1 }, { "geno1": 2, "geno2": 2 }
    ],
    "distance": 1000
  }
}
```

**Response:**
```json
{ "variant1": "V1", "variant2": "V2", "dPrime": 0.0, "rSquared": 0.0, "distance": 1000 }
```

## Performance

- **Time Complexity:** O(n) in the number of individuals
- **Space Complexity:** O(1)

## References

- Hill, W. G., Robertson, A. (1968). Linkage disequilibrium in finite populations. *Theor. Appl. Genet.* 38:226–231.
- Lewontin, R. C. (1964). The interaction of selection and linkage. *Genetics* 49:49–67.

## See Also

- [haplotype_blocks](haplotype_blocks.md)
