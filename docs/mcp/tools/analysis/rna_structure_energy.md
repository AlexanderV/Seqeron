# rna_structure_energy

Turner 2004 free energy of a given RNA secondary structure (RNAeval counterpart).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_structure_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateFreeEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Evaluates the free energy ΔG°37 (kcal/mol) of a pseudoknot-free dot-bracket structure loop by loop with
exactly the energy model minimised by `minimum_free_energy`: exterior helix ends with d2 terminal mismatches /
dangles and AU/GU penalties, hairpins (special loops, mismatch bonuses, all-C), stacks (incl. the GGUC/CUGG
three-stack), bulges, interior loops (int11/int21/int22 and the generic rules) and linear d2 multiloops
(9.3 − 0.9 per helix). It is the counterpart of ViennaRNA 2.7.2 `RNAeval -d2` / `eval_structure`; the energy of
the MFE structure equals the MFE.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) (`CalculateFreeEnergy`)
- Algorithm: [RNA_Free_Energy.md](../../../algorithms/RnaStructure/RNA_Free_Energy.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (T read as U, case-insensitive) |
| `dotBracket` | string | Yes | Structure of the same length using only `(`, `)` and `.` |
| `temperature` | number | No | Absolute temperature in Kelvin (default 310.15 = 37 °C); every Turner 2004 parameter is rescaled with its enthalpy, dG(T) = dH - (dH - dG37)*T/310.15 (ViennaRNA md.temperature); must be positive |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Free energy of the structure (kcal/mol, 2 dp) |

## Errors

Empty sequence; null structure; a character other than `(`, `)`, `.`; unbalanced or wrong-length structure;
a pair of bases that cannot pair; a hairpin loop of fewer than 3 unpaired nucleotides (prohibited by the Turner 2004
rules; ViennaRNA `eval_structure` reports INF).

## Examples

### Example 1: hairpin

```json
{ "tool": "rna_structure_energy", "arguments": { "rnaSequence": "GGGAAACCC", "dotBracket": "(((...)))" } }
```

**Response:** `{ "energy": -1.12 }` (ViennaRNA `eval_structure`, NNDB-precision parameters: −1.12).

### Example 2: three-way multiloop

```json
{ "tool": "rna_structure_energy", "arguments": { "rnaSequence": "GAGGGAAAACCCAGGGAAAACCCAC", "dotBracket": "(.(((....))).(((....))).)" } }
```

**Response:** `{ "energy": -1.54 }` (ViennaRNA: −1.54).

## Performance

- **Time Complexity:** O(n) loops, O(n) per loop scan.
- **Space Complexity:** O(n).

## See Also

- [minimum_free_energy](minimum_free_energy.md)
- [multibranch_loop_energy](multibranch_loop_energy.md)
- [terminal_mismatch_energy](terminal_mismatch_energy.md)
