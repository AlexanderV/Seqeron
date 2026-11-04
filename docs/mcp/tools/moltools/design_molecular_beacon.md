# design_molecular_beacon

Design a hairpin molecular-beacon probe for real-time detection.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_molecular_beacon` |
| **Method ID** | `ProbeDesigner.DesignMolecularBeacon` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Chooses the highest-scoring `probe_length`-bp loop in the target (GC 40–60%, loop Tm in the window, no homopolymer run > 4 preferred) and flanks it with GC-rich complementary stems to form a hairpin. The 5′ stem is `⌊stem_length/2⌋` G's followed by the remaining C's; the 3′ stem is its reverse complement. The final probe is `stem5 + loop + stem3`; `Tm` is the loop (probe–target) Tm — Primer3 `seqtm` at the Primer3 probe conditions (50 nM, 50 mM monovalent, no Mg²⁺/dNTP) — and `Start`/`End` mark the loop in the target. The warnings carry the ntthal stem-loop (hairpin) Tm of the whole beacon (beacons ≤ `max_align_length` nt, default 60 = Primer3's THAL_MAX_ALIGN; opt-in up to 10 000 = thal.c compiled with `-DTHAL_MAX_ALIGN`). With `detection_temperature` = T the loop Tm window is [T + 7, T + 10] °C and the stem-loop Tm is checked against T + 7 °C (molecular-beacon rules of Tyagi & Kramer: probe and stem Tm 7–10 °C above the detection temperature); without it the window is 55–65 °C. Returns `probe = null` when the target is shorter than `probe_length`.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1964](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1964)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_sequence` | string | Yes | Target DNA sequence, non-empty. |
| `probe_length` | integer | No | Loop length in bp (> 0, default 25). |
| `stem_length` | integer | No | Stem length in bp (> 0, default 5). |
| `detection_temperature` | number | No | Detection (annealing) temperature in °C for the 7–10 °C rules (default: none). |
| `max_align_length` | integer | No | THAL_MAX_ALIGN of the ntthal stem-loop Tm (60–10 000, default 60 = Primer3): longest beacon that gets it. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probe` | object \| null | `Probe` (`Type = MolecularBeacon`) or null if the target is too short. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Target sequence cannot be null or empty |
| 1002 | Probe (loop) length must be positive |
| 1003 | Stem length must be positive |

## Examples

### Example 1: 48-nt target, 20 bp loop, 5 bp stem

`stem5 = "GGCCC"`, `stem3 = revComp("GGCCC") = "GGGCC"`; beacon = `GGCCC + target[0..19] + GGGCC` = 30 bp, `Type = MolecularBeacon`, loop `Start = 0`, `End = 19`.

**Input:** `{ "target_sequence": "ACGTACGT…", "probe_length": 20, "stem_length": 5 }`

**Response (abridged):** `{ "probe": { "sequence": "GGCCC…GGGCC", "start": 0, "end": 19, "type": "MolecularBeacon" } }`

### Example 2: Too-short target

`design_molecular_beacon("ACGT", 20)` returns `{ "probe": null }`.

## See Also

- [design_probes](design_probes.md), [validate_probe](validate_probe.md)
