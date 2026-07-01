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

Chooses the highest-scoring `probe_length`-bp loop in the target (GC 40–60%, Tm 55–65 °C, no long homopolymer preferred) and flanks it with GC-rich complementary stems to form a hairpin. The 5′ stem is `⌊stem_length/2⌋` G's followed by the remaining C's; the 3′ stem is its reverse complement. The final probe is `stem5 + loop + stem3`; `Tm` is the loop Tm and `Start`/`End` mark the loop in the target. Returns `probe = null` when the target is shorter than `probe_length`.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L785](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L785)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_sequence` | string | Yes | Target DNA sequence, non-empty. |
| `probe_length` | integer | No | Loop length in bp (> 0, default 25). |
| `stem_length` | integer | No | Stem length in bp (> 0, default 5). |

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
