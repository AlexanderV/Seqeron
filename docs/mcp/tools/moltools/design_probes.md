# design_probes

Design ranked hybridization probes for a target sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_probes` |
| **Method ID** | `ProbeDesigner.DesignProbes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans the target for every candidate of admissible length (`parameters.MinLength..MaxLength`), scores each on GC%, Tm, homopolymers, self-complementarity and structure heuristics, and returns up to `max_probes` probes sorted by score (descending). Use a `ProbeParameters` preset (`Microarray` default, `FISH`, `NorthernBlot`, `qPCR`, `SouthernBlot`) or custom values. A target shorter than the minimum probe length returns an empty list.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L493](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L493)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_sequence` | string | Yes | Target DNA sequence, non-empty. |
| `parameters` | object | No | Optional `ProbeParameters` (defaults to Microarray). |
| `max_probes` | integer | No | Maximum probes to return (> 0, default 10). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | array | Up to `max_probes` `Probe` records (sequence, start, end, Tm, GC, score, type, warnings), score-descending. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Target sequence cannot be null or empty |
| 1002 | Maximum probes must be positive |

## Examples

### Example 1: 81-nt ACGT repeat, top 3 probes

**Input:** `{ "target_sequence": "ACGTACGT…ACGT", "max_probes": 3 }`

Returns exactly 3 Microarray probes (length 50–60), each a substring of the target, ordered by descending score.

### Example 2: Too-short target

`design_probes("ACGTACGTACGT")` returns `{ "probes": [] }` (shorter than the 50-nt Microarray minimum).

## See Also

- [design_tiling_probes](design_tiling_probes.md), [design_antisense_probes](design_antisense_probes.md), [validate_probe](validate_probe.md)
