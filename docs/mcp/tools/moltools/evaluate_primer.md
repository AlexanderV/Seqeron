# evaluate_primer

Evaluate a single primer against quality criteria.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `evaluate_primer` |
| **Method ID** | `PrimerDesigner.EvaluatePrimer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scores a single primer and returns a candidate record: length, GC%, Tm (Primer3-default SantaLucia 1998 NN Tm: 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM), longest homopolymer, hairpin potential, 3′-end ΔG°37 stability, a list of quality issues (against the supplied or default `PrimerParameters`), an overall validity flag, an informational numeric score and the Primer3 per-primer `penalty` (|Tm − OptimalTm| + |length − OptimalLength|; lower is better). `position` and `is_forward` are informational and echoed back.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L120](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L120)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Primer sequence (non-empty). |
| `position` | integer | Yes | 0-based location (informational). |
| `is_forward` | boolean | Yes | Forward/reverse flag. |
| `parameters` | object | No | Optional design parameters. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` / `gcContent` / `meltingTemperature` | number | Basic metrics. |
| `homopolymerLength` / `hasHairpin` / `stability3Prime` | mixed | Structural metrics. |
| `isValid` / `issues` / `score` / `penalty` | mixed | QC verdict; `penalty` is the Primer3 ranking penalty. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATCGATCGATCGATCGATCG`

20-mer, 50% GC, Primer3-default Tm = 57.36 °C (primer3-py `calc_tm`) → `57.4`, penalty 2.637, homopolymer length 1.

## See Also

- [design_primers](design_primers.md), [primer_melting_temperature](primer_melting_temperature.md), [three_prime_stability](three_prime_stability.md)
