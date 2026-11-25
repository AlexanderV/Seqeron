# rna_pseudoknot_structure

MFE RNA structure that may contain canonical H-type pseudoknots (pknotsRG class).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_pseudoknot_structure` |
| **Method ID** | `RnaSecondaryStructure.PredictStructurePseudoknotRecursive` (`allowNestedKnots = true`), `RnaSecondaryStructure.PredictStructurePseudoknot` (`false`) |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Exact optimum of the Reeder & Giegerich (2004) canonical simple recursive pseudoknot grammar (class csr-PK) over the
shared Turner 2004 / ViennaRNA d2 model: equal-length, bulge-free, maximal crossing helices a·a' / b·b' with three loops
u, v, w, pseudoknot initiation 9.0, 0.3 kcal/mol per unpaired knot-loop nucleotide, 6.0 for a knot used as a multiloop or
knot-loop component. `allowNestedKnots = false` allows at most one knot, in the exterior loop, with pseudoknot-free
loops; `true` (default) allows several knots, knots in multiloops and knots in knot loops. The result is never above the
pseudoknot-free MFE (`rna_mfe_structure`) and equals it when no knot pays off. O(n^4) time, O(n^2) space.
Details, deviations from pknotsRG (Turner 2004 instead of Turner 1999 parameters; no relaxed 1-nt-bulge helices) and the
oracle: [Pseudoknot_Prediction.md](../../../algorithms/RnaStructure/Pseudoknot_Prediction.md),
[Pseudoknot_Prediction_Recursive.md](../../../algorithms/RnaStructure/Pseudoknot_Prediction_Recursive.md).

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence; T read as U, case-insensitive; other characters never pair |
| `allowNestedKnots` | boolean | No | `true` (default) full recursive grammar; `false` at most one knot in the exterior loop |
| `minLoopSize` | integer | No | Minimum hairpin loop size (default 3; values < 3 behave as 3) |
| `temperature` | number | No | Kelvin (default 310.15); Turner terms rescaled with their enthalpies; pseudoknot constants are temperature independent |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Upper-cased sequence (T → U) |
| `dotBracket` | string | Structure; crossing helices on further bracket families (`[]`, `{}`, ...) |
| `basePairs` | array | Base pairs (0-based, sorted) |
| `freeEnergy` | number | ΔG (kcal/mol) |
| `hasPseudoknot` | boolean | The pair set contains a crossing pair |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

```json
{ "tool": "rna_pseudoknot_structure", "arguments": { "rnaSequence": "GGGGAACCCCAACCCCAAGGGG" } }
```
→ `((((..[[[[..))))..]]]]`, −13.56, `hasPseudoknot: true` (plain MFE −7.84).

## See Also

- [rna_mfe_structure](rna_mfe_structure.md) — pseudoknot-free optimum
- [detect_pseudoknots](detect_pseudoknots.md), [pseudoknot_orders](pseudoknot_orders.md)
