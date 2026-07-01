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

Scores a single primer and returns a candidate record: length, GC%, Tm (Wallace/Marmur–Doty), longest homopolymer, hairpin potential, 3′-end ΔG°37 stability, a list of quality issues (against the supplied or default `PrimerParameters`), an overall validity flag, and a numeric score. `position` and `is_forward` are informational and echoed back.

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
| `isValid` / `issues` / `score` | mixed | QC verdict. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATCGATCGATCGATCGATCG`

20-mer, 50% GC, Marmur–Doty Tm = `64.9 + 41·(10−16.4)/20 = 51.8` °C, homopolymer length 1.

## See Also

- [design_primers](design_primers.md), [primer_melting_temperature](primer_melting_temperature.md), [three_prime_stability](three_prime_stability.md)
