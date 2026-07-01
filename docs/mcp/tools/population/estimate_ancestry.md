# estimate_ancestry

Estimate per-individual ancestry proportions against fixed reference populations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `estimate_ancestry` |
| **Method ID** | `PopulationGeneticsAnalyzer.EstimateAncestry` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Supervised (projection) ADMIXTURE: for each individual the ancestry vector q (one fraction per reference population) is found by maximum likelihood with the reference allele-1 frequencies held fixed, using the FRAPPE expectation-maximization update on the binomial admixture log-likelihood (Alexander, Novembre & Lange 2009, Eqs. 2/4/5).

- Genotype `Genotypes[j]` is the number of copies of allele 1 at SNP j (0, 1, or 2); any other value is treated as missing and that SNP is skipped.
- Each individual's genotype length must equal the reference-panel SNP count; mismatched individuals are skipped.
- Proportions are keyed by reference-population id and sum to 1. With no informative SNP (or `maxIterations = 0`) the uniform prior 1/K is returned.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L1264](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L1264)
- Algorithm doc: [Ancestry_Estimation.md](../../../algorithms/Population_Genetics/Ancestry_Estimation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `individuals` | object[] | Yes | Individuals with `individualId` and `genotypes` (0/1/2 per SNP) |
| `referencePops` | object[] | Yes | Reference populations with `populationId` and `alleleFrequencies` (allele-1 freq per SNP) |
| `maxIterations` | integer | No | Maximum EM iterations per individual (default 100) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | One entry per valid individual |
| `items[].individualId` | string | Individual identifier |
| `items[].proportions` | object | Map of reference-population id → ancestry fraction (sums to 1) |

## Errors

None (empty individuals or empty reference panels return an empty result).

## Examples

### Example 1: One EM iteration on a symmetric panel

**User Prompt:**
> Run one iteration of ancestry estimation for genotype (2, 0) against panels A=(0.8,0.2), B=(0.2,0.8).

**Expected Tool Call:**
```json
{
  "tool": "estimate_ancestry",
  "arguments": {
    "individuals": [{ "individualId": "IND1", "genotypes": [2, 0] }],
    "referencePops": [
      { "populationId": "A", "alleleFrequencies": [0.8, 0.2] },
      { "populationId": "B", "alleleFrequencies": [0.2, 0.8] }
    ],
    "maxIterations": 1
  }
}
```

**Response:**
```json
{ "items": [{ "individualId": "IND1", "proportions": { "A": 0.8, "B": 0.2 } }] }
```

Each SNP mix = 0.5; Eq. 4 gives q_A = (1.6 + 1.6)/(2·2) = 0.8, q_B = 0.2.

### Example 2: Converged diagnostic individual

**User Prompt:**
> Estimate ancestry to convergence for the same input.

**Expected Tool Call:**
```json
{
  "tool": "estimate_ancestry",
  "arguments": {
    "individuals": [{ "individualId": "IND1", "genotypes": [2, 0] }],
    "referencePops": [
      { "populationId": "A", "alleleFrequencies": [0.8, 0.2] },
      { "populationId": "B", "alleleFrequencies": [0.2, 0.8] }
    ],
    "maxIterations": 1000
  }
}
```

**Response:**
```json
{ "items": [{ "individualId": "IND1", "proportions": { "A": 1.0, "B": 0.0 } }] }
```

The EM ascends to the MLE: this genotype is fully population-A ancestry.

## Performance

- **Time Complexity:** O(individuals · iterations · SNPs · K)
- **Space Complexity:** O(K) per individual

## References

- Alexander, D. H., Novembre, J., Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 19(9):1655–1664.

## See Also

- [f_statistics](f_statistics.md), [pairwise_fst](pairwise_fst.md)
