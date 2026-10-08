# design_tiling_probes

Tile a target with fixed-length, overlapping probes for full coverage.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_tiling_probes` |
| **Method ID** | `ProbeDesigner.DesignTilingProbes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Walks the target in steps of `probe_length − overlap`, emitting a `probe_length`-bp probe at each start `0, step, 2·step, …` up to `target.Length − probe_length`; when that grid stops short of the 3′ end (`(target.Length − probe_length) mod step ≠ 0`) one more probe is anchored at the end (start `target.Length − probe_length`), as CATCH (Metsky et al. 2019) does, so every base is covered (B07 audit round 7, A7-2: previously the tail after the last grid probe was left uncovered). The target must be at least `probe_length` long. Sub-optimal probes are still emitted (with a `"Suboptimal probe, included for coverage"` warning) so coverage is preserved. Returns the probe list plus the number of covered positions, the mean probe Tm, and the Tm range (max − min).

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1319](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1319)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_sequence` | string | Yes | Target DNA sequence, non-empty. |
| `probe_length` | integer | No | Probe length in bp (> 0, default 60). |
| `overlap` | integer | No | Overlap in bp (`0 ≤ overlap < probe_length`, default 20). |
| `parameters` | object | No | Optional `ProbeParameters` (Tm/GC bounds etc.). |
| `thermodynamic_screen_max_length` | integer | No | THAL_MAX_ALIGN of the ntthal self-structure screen (`ProbeParameters.ThermodynamicScreenMaxLength`): A/C/G/T probes up to this length get the ntthal screen. Default null = the parameters' value (60 = Primer3); opt-in 61–10 000 screens longer probes with the unchanged ntthal recursions (= thal.c compiled with `-DTHAL_MAX_ALIGN`; cost O(n²) per probe). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | array | Tiling `Probe` records (all `Type = Tiling`). |
| `coverage` | integer | Number of target positions covered (= target length). |
| `meanTm` | number | Mean probe Tm. |
| `tmRange` | number | max(Tm) − min(Tm). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Target sequence cannot be null or empty |
| 1002 | Probe length must be positive |
| 1003 | Overlap must be non-negative and less than the probe length |
| 1004 | Probe length must not exceed the target length |

## Examples

### Example 1: 208-nt target, 50 bp probes, 10 bp overlap

Target `A×100 + GCGCGCGC + T×100`; step = 40 → grid starts `{0, 40, 80, 120}` (covering 0–169) plus the end-anchored probe at 208 − 50 = 158 → 5 probes covering all 208 positions.

**Input:** `{ "probe_length": 50, "overlap": 10 }`

**Response (abridged):** `{ "coverage": 208, "probes": [ { "start": 0 }, { "start": 40 }, { "start": 80 }, { "start": 120 }, { "start": 158 } ] }`

### Example 2: Invalid overlap

`design_tiling_probes("…", 50, 50)` (overlap ≥ length) throws `ArgumentException`; so does a target shorter than `probe_length` (e.g. 30 nt with the default 60).

## See Also

- [design_probes](design_probes.md), [validate_probe](validate_probe.md)
