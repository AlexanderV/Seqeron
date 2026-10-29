# calculate_on_target_doench2014

Doench et al. 2014 "Rule Set 1" on-target efficacy score for an SpCas9 guide (0..100, unrounded).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `calculate_on_target_doench2014` |
| **Method ID** | `CrisprDesigner.CalculateOnTargetDoench2014` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Doench et al. 2014 "Rule Set 1" on-target efficacy score for an SpCas9 guide, returned on a 0..100 scale (higher = predicted more active). This is the published logistic linear model: intercept `0.59763615` + a GC-count term over the 20-nt protospacer using the paper's **two** coefficients (`-0.2026259` for GC <= 10, `-0.1665878` for GC > 10, applied as `abs(10 - GC) x coefficient`) + 70 position-specific single- and di-nucleotide weights, through a logistic sigmoid.

Source: Doench et al. 2014, Nat Biotechnol 32:1262 (PMID 25184501); coefficients taken verbatim from the reference implementation `doenchScore.py` shipped in CRISPOR (Haeussler et al. 2016, Genome Biol 17:148).

Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be NGG) + 3 nt downstream, A/C/G/T only (case-insensitive); the NGG requirement is this library's explicit guard (the model's own feature table presumes NGG).

Returned **unrounded**: CRISPOR's own table reports `int(100 x sigmoid)`, i.e. the floor of this value. Call to rank guides by predicted cutting efficiency when the 30-nt context is known. This is **not** Rule Set 2 / Azimuth — for that, use [calculate_on_target_rule_set2](calculate_on_target_rule_set2.md).

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1284](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1284)
- Algorithm doc: [Guide_RNA_Design.md](../../../algorithms/MolTools/Guide_RNA_Design.md) §2.3 "Published on-target efficacy models"

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `context_30mer` | string | Yes | 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt NGG PAM + 3 nt downstream. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Predicted on-target activity, 0..100, unrounded (CRISPOR's column is its floor). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Context 30-mer cannot be null or empty. |
| 1002 | Doench 2014 requires a 30-nt context (4 upstream + 20 protospacer + 3 PAM + 3 downstream); got {length}. |
| 1003 | Doench 2014 context must contain only A/C/G/T; found '{char}'. |
| 1004 | Doench 2014 expects an SpCas9 NGG PAM at offsets 25-26 of the 30-nt context. |

## Examples

### Example 1: first worked example of the reference

`context_30mer: "TATAGCTGCGATCTGAGGTAGGGAGGGACC"` → `score: 71.30893469547244` (CRISPOR column: `71`).

### Example 2: second worked example of the reference

`context_30mer: "TCCGCACCTGTCACGGTCGGGGCTTGGCGC"` → `score: 1.8983843193336125` (CRISPOR column: `1`).

### Example 3: no NGG PAM

`context_30mer: "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"` → error 1004.

## Performance

- **Time Complexity:** O(1) — fixed 30-nt input, 70 weights.
- **Space Complexity:** O(1).

## See Also

- [calculate_on_target_rule_set2](calculate_on_target_rule_set2.md), [design_guide_rnas](design_guide_rnas.md), [evaluate_guide_rna](evaluate_guide_rna.md)
