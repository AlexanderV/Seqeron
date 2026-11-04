# design_antisense_probes

Design antisense probes against an mRNA-sense sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_antisense_probes` |
| **Method ID** | `ProbeDesigner.DesignAntisenseProbes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reverse-complements the supplied mRNA-sense sequence and runs the probe designer on it, returning up to `max_probes` top-scoring probes each tagged `Type = Antisense`. Each probe sequence is therefore a substring of the reverse complement of the input mRNA. Uses the same `ProbeParameters` presets as `design_probes` (default Microarray).

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1401](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1401)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mrna_sequence` | string | Yes | mRNA-sense sequence (reverse-complemented internally), non-empty. |
| `parameters` | object | No | Optional `ProbeParameters` (defaults to Microarray). |
| `max_probes` | integer | No | Maximum probes to return (> 0, default 5). |
| `thermodynamic_screen_max_length` | integer | No | THAL_MAX_ALIGN of the ntthal self-structure screen (`ProbeParameters.ThermodynamicScreenMaxLength`): A/C/G/T probes up to this length get the ntthal screen. Default null = the parameters' value (60 = Primer3); opt-in 61–10 000 screens longer probes with the unchanged ntthal recursions (= thal.c compiled with `-DTHAL_MAX_ALIGN`; cost O(n²) per probe). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | array | `Probe` records with `Type = Antisense`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1002 | Maximum probes must be positive |

## Examples

### Example 1: 81-nt mRNA repeat

**Input:** `{ "mrna_sequence": "ACGTACGT…ACGT", "max_probes": 5 }`

Returns 5 antisense probes; each `probe.Sequence` equals a substring of `reverseComplement(mrna)`.

### Example 2: Empty input

`design_antisense_probes("")` throws `ArgumentException`.

## See Also

- [design_probes](design_probes.md), [design_tiling_probes](design_tiling_probes.md)
