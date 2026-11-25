# rna_region_accessibility

Probability that a contiguous RNA region is entirely unpaired in the equilibrium ensemble (structural accessibility), with its opening energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_region_accessibility` |
| **Method ID** | `RnaSecondaryStructure.CalculateRegionUnpairedProbability` (+ `CalculateRegionOpeningEnergy`) |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns `P(region unpaired) = Z_open / Z` of the McCaskill partition function on the Turner 2004 model (the engine of
`rna_partition_function`), where `Z_open` sums the Boltzmann weights of the structures in which no base of the region
pairs, and the opening energy `-RT ln P` (kcal/mol; Kertesz et al. 2007). The ensemble is the whole sequence: the value
equals ViennaRNA `RNAplfold -u` for window `W >= n` and span `L >= n` (`RNA.pfl_fold_up(seq, u, n, n)`) and the exact
hard-constraint ratio `exp(-(G_constrained - G_ensemble)/RT)`. It is not the sliding-window local average of RNAplfold with
`W < n`, and it is not the product of the per-base unpaired probabilities returned by `rna_partition_function`.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.PartitionFunction.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.PartitionFunction.cs)
- Algorithm: [Turner_McCaskill_Partition_Function.md](../../../algorithms/RnaStructure/Turner_McCaskill_Partition_Function.md)
- Evidence: [RNA-ACCESS-001-Evidence.md](../../../Evidence/RNA-ACCESS-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (T read as U, case-insensitive; N / IUPAC never pair) |
| `windowEnd` | integer | Yes | 0-based inclusive last position of the region |
| `windowLength` | integer | Yes | Region length u >= 1 (TargetScan uses 14); the region must lie inside the sequence |
| `minLoopSize` | integer | No | Minimum hairpin loop size (default 3; values < 3 clamped) |
| `temperature` | number | No | Kelvin (default 310.15 = 37 °C); must be positive |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `windowStart`, `windowEnd`, `windowLength` | integer | The queried region (0-based inclusive) |
| `unpairedProbability` | number | Z_open / Z in [0, 1] |
| `openingEnergy` | number or null | -RT ln P (kcal/mol); null when P underflows to 0 |
| `temperature` | number | Kelvin |

## Errors

Empty sequence; `windowLength < 1`; region outside the sequence (including `windowEnd` negative or >= n); temperature <= 0 or non-finite.

## Examples

```json
{ "tool": "rna_region_accessibility",
  "arguments": { "rnaSequence": "GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", "windowEnd": 34, "windowLength": 14 } }
```

**Response (abridged):** `unpairedProbability 0.0305` (ViennaRNA constrained ensemble 0.030532), `openingEnergy 2.15`.

## Performance

- **Time Complexity:** O(n^3) (two partition-function fills). **Space Complexity:** O(n^2).

## See Also

- [rna_partition_function](rna_partition_function.md)
- [rna_mfe_structure](rna_mfe_structure.md)
