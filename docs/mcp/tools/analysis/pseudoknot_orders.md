# pseudoknot_orders

Assign pseudoknot orders to RNA base pairs and render dot-bracket-letter notation.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `pseudoknot_orders` |
| **Method ID** | `RnaSecondaryStructure.AssignPseudoknotOrders` (+ `ToDotBracketByPseudoknotOrder`, `ToDotBracket`) |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Order 0 is a **maximum-cardinality nested** (pseudoknot-free) subset of the pairs; the removed pairs are
processed again for order 1, and so on (biotite `structure.pseudoknots`, the Smit et al. 2008 optimal
removal applied iteratively). Pairs of equal order never cross. When several maximum subsets exist the one
keeping the 5'-most opening pair is chosen. `orderDotBracket` writes order *k* with the *k*-th family of
`() [] {} <> Aa Bb …` (Antczak et al. 2018 DBL notation). `dotBracket` is the ViennaRNA 2.7.2
`vrna_db_from_ptable` rendering (5'→3' first-fit passes over the same families) — it can differ from the
order string, e.g. when the larger helix opens 3' of a smaller crossing one.

The number of orders is an upper bound on the minimum decomposition count (the chromatic number of the
crossing graph, NP-complete to compute in general).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L2199](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L2199)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `basePairs` | object[] | Yes | Base pairs: `{ position1, position2, base1, base2, type }` (0-based, either orientation) |
| `length` | integer | Yes | Structure length; every position must be `< length` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `orders` | int[] | Order of each input pair (index-aligned) |
| `orderDotBracket` | string | Dot-bracket-letter string by order |
| `dotBracket` | string | ViennaRNA `vrna_db_from_ptable` rendering |

## Errors

| Code | Message |
|------|---------|
| 1002 | Unknown base pair type: … |
| 1001 | Self pair, position outside `[0, length)`, a position in two pairs, or more than 30 families needed (ViennaRNA would silently drop the excess pairs) |

## Examples

### Example 1: biotite docstring example "ABCbac"

**Input:** pairs `(0,4)`, `(1,3)`, `(2,5)`, `length` 6 →
**Response:** `{ "orders": [0,0,1], "orderDotBracket": "(([))]", "dotBracket": "(([))]" }`

### Example 2: larger 3' helix

**Input:** pairs `(0,6)`, `(2,10)`, `(3,9)`, `(4,8)`, `length` 12 →
`orders` `[1,0,0,0]`, `orderDotBracket` `[.(((.].))).`, `dotBracket` `(.[[[.).]]].` (= `RNA.db_from_ptable`).

## Performance

- **Time:** O(K·m²) for m = 2·pairs and K orders; **Space:** O(m²) (interval DP per layer). `dotBracket`: O(F·n).

## See Also

- [detect_pseudoknots](detect_pseudoknots.md)
- [parse_dot_bracket](parse_dot_bracket.md)
