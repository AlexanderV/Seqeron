# rna_mfe_structure

Minimum-free-energy RNA secondary structure (RNAfold counterpart).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_mfe_structure` |
| **Method ID** | `RnaSecondaryStructure.PredictStructureMfe` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Folds an RNA sequence into its minimum-free-energy pseudoknot-free secondary structure by Zuker–Stiegler
dynamic programming and traceback with the Turner 2004 nearest-neighbour parameters and the ViennaRNA d2 model
(`RNAfold` default: terminal mismatches / dangles on every helix end, linear multiloop 9.3 + 0·unpaired − 0.9·branch,
interior loops ≤ 30 nt, lonely pairs allowed). Returns the dot-bracket, the base pairs and the MFE, which equals
`minimum_free_energy` and `rna_structure_energy` of the returned structure. When several structures share the
optimum the traceback returns one of them deterministically (it can differ from RNAfold's choice; the energy is the
same). Unlike `predict_rna_structure` (assembly of isolated hairpin stem-loops), this is the global optimum including
multiloops, bulges and interior loops.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) (`PredictStructureMfe`, `CalculateMfeStructure`)
- Algorithm: [Minimum_Free_Energy.md](../../../algorithms/RnaStructure/Minimum_Free_Energy.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence; T read as U, case-insensitive; other characters (N, IUPAC, gaps) never pair |
| `minLoopSize` | integer | No | Minimum hairpin loop size (default 3; values < 3 behave as 3, as ViennaRNA `md.min_loop_size`) |
| `temperature` | number | No | Absolute temperature in Kelvin (default 310.15 = 37 °C); every Turner 2004 parameter is rescaled with its enthalpy, dG(T) = dH - (dH - dG37)*T/310.15 (ViennaRNA md.temperature); must be positive |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Upper-cased sequence (T → U) |
| `dotBracket` | string | MFE structure |
| `basePairs` | array | Base pairs (0-based) |
| `minimumFreeEnergy` | number | ΔG°37 of the structure (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Hairpin

**User Prompt:**
> Fold "GGGAAAACCC".

**Expected Tool Call:**
```json
{ "tool": "rna_mfe_structure", "arguments": { "rnaSequence": "GGGAAAACCC" } }
```

**Response:**
```json
{ "sequence": "GGGAAAACCC", "dotBracket": "(((....)))", "minimumFreeEnergy": -2.02 }
```
ViennaRNA 2.7.2 `RNA.fold` (NNDB-precision Turner 2004 file): `(((....)))` −2.02.

### Example 2: Two hairpins

```json
{ "tool": "rna_mfe_structure", "arguments": { "rnaSequence": "GAGGGAAAACCCAGGGAAAACCCAC" } }
```
→ `..(((....))).(((....)))..`, −7.04 (= `RNA.fold`).

## Performance

- **Time Complexity:** O(n³) (multiloop split) + O(n²·30²) (interior loops); ~0.8 s at n = 500, ~3.5 s at n = 1000 (Release).
- **Space Complexity:** O(n²).

## See Also

- [minimum_free_energy](minimum_free_energy.md) — the MFE value only
- [rna_structure_energy](rna_structure_energy.md) — energy of a given structure
- [predict_rna_structure](predict_rna_structure.md) — stem-loop assembly
