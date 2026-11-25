# rna_partition_function

McCaskill equilibrium partition function of an RNA on the Turner 2004 model (ViennaRNA `RNAfold -p` counterpart).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_partition_function` |
| **Method ID** | `RnaSecondaryStructure.CalculatePartitionFunction` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sums the Boltzmann weights of all pseudoknot-free secondary structures with the same Turner 2004 loop model
(d2 dangles) as `minimum_free_energy` and returns the partition function `Z`, the ensemble free energy `-RT ln Z`
(ViennaRNA `pf()`), the base-pair probabilities (`bpp()`), the per-base unpaired probabilities, the MFE structure and its
Boltzmann probability (`pr_structure`), the ensemble diversity (`mean_bp_distance`) and the centroid structure (`centroid`).
`temperature` rescales every Turner parameter with its enthalpy (`dG(T) = dH - (dH - dG37)*T/310.15`, ViennaRNA
`md.temperature`), RT = 1.98717e-3·T kcal/mol. The matrices are scaled per nucleotide (`pf_scale`), so sequences of
thousands of nucleotides neither overflow nor underflow; `partitionFunction` is `null` when Z exceeds the double range
(use `logPartitionFunction`).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.PartitionFunction.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.PartitionFunction.cs)
- Algorithm: [RNA_Partition_Function.md](../../../algorithms/RnaStructure/RNA_Partition_Function.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (T read as U, case-insensitive; N / IUPAC never pair) |
| `minLoopSize` | integer | No | Minimum hairpin loop size (default 3; values < 3 clamped) |
| `temperature` | number | No | Kelvin (default 310.15 = 37 °C); must be positive |
| `minProbability` | number | No | List only pairs with P >= this value (default 0.001; 0 lists every pair with P > 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `partitionFunction` | number or null | Z (null beyond the double range) |
| `logPartitionFunction` | number | ln Z |
| `ensembleFreeEnergy` | number | -RT ln Z (kcal/mol) |
| `basePairProbabilities` | array | `{position1, position2, probability}` (0-based, position1 < position2) |
| `unpairedProbabilities` | array | P(unpaired) per position |
| `mfeStructure`, `mfeFreeEnergy`, `mfeProbability` | string, number, number | MFE structure, its energy and Boltzmann probability |
| `ensembleDiversity` | number | mean base-pair distance |
| `centroidStructure` | string | pairs with P > 0.5 |
| `temperature` | number | Kelvin |

## Errors

Empty sequence; temperature <= 0 or non-finite; `minProbability` outside [0,1].

## Examples

```json
{ "tool": "rna_partition_function", "arguments": { "rnaSequence": "GGGAAAACCC" } }
```

**Response (abridged):** `mfeStructure "(((....)))"`, `mfeFreeEnergy -2.02`, `ensembleFreeEnergy -2.1535` (ViennaRNA pf -2.1529),
`mfeProbability 0.805`, `ensembleDiversity 1.32`.

## Performance

- **Time Complexity:** O(n^3). **Space Complexity:** O(n^2) (Release: n = 500 1.3 s, n = 1000 6.5 s, n = 2000 52 s).

## See Also

- [rna_mfe_structure](rna_mfe_structure.md)
- [minimum_free_energy](minimum_free_energy.md)
- [rna_structure_energy](rna_structure_energy.md)
