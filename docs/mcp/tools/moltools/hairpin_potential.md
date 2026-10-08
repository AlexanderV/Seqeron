# hairpin_potential

Detect whether a sequence can fold into a hairpin.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `hairpin_potential` |
| **Method ID** | `PrimerDesigner.HasHairpinPotential` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns true if the sequence contains an exact Watson–Crick self-complementary stem of at least `min_stem_length` bp separated by a loop of at least `min_loop_length` nt (no G·T wobble, mismatches or energies). A short-sequence O(n²) scan is used below 100 bp; a suffix-tree scan is used at or above 100 bp. Sequences shorter than `2·min_stem_length + min_loop_length` cannot form a hairpin and return false.

**Library screen (sequence-only heuristic, not Primer3's).** The default minimum loop of 3 nt is Primer3's
hairpin minimum (`thal.c`: `static const int min_hrpn_loop = 3;`); the default minimum stem of 4 bp is an
unsourced library threshold — no published definition of a "≥ 4-bp stem + ≥ 3-nt loop" rule was found
(audit round 3, A3-8). The screen can disagree with Primer3's thermodynamic hairpin screen (primer3-py 2.3.1
`calc_hairpin`, 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM):

| Sequence | `hasHairpin` | Primer3 hairpin |
|----------|--------------|-----------------|
| `AAAACCCTTTT` | true | no structure found |
| `CAGTAAAACCCTTTTGCAGC` | true | Tm 37.65 °C (< `PRIMER_MAX_HAIRPIN_TH` 47 °C, passes) |

For Primer3's sourced hairpin screen use [evaluate_primer](evaluate_primer.md) (`hairpinTh` = ntthal hairpin
Tm, rejected above 47 °C) or the C# API `PrimerDesigner.CalculateHairpinThermodynamicsNtthal` /
`CalculatePrimer3OligoStructure`.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L1786](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L1786)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (non-empty). |
| `min_stem_length` | integer | No | Minimum stem length (default 4, positive). |
| `min_loop_length` | integer | No | Minimum loop length (default 3, non-negative). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hasHairpin` | boolean | True if a hairpin can form. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Minimum stem length must be positive |
| 1003 | Minimum loop length cannot be negative |

## Examples

### Example 1: `GGGGAAACCCC` → `true` (stem `GGGG` / loop `AAA` / stem `CCCC`).

### Example 2: `AAAAAAAAAAA` → `false` (no complementary stem).

## See Also

- [three_prime_stability](three_prime_stability.md), [evaluate_primer](evaluate_primer.md)
