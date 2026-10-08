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

Scans the target for every candidate of admissible length (`parameters.MinLength..MaxLength`), scores each with an additive penalty score on GC%, Tm (Primer3 `seqtm` at the parameters' conditions; the default Microarray preset uses OligoArray 2.0's conditions — 1 M Na⁺, 1 µM oligo, nearest-neighbour Tm over the whole ≤ 60-nt probe, Tm window 82–90 °C; the other presets use the Primer3 probe conditions 50 nM / 50 mM monovalent / no Mg²⁺ / no dNTP), homopolymers, self-structure (for ACGT probes of ≤ `thermodynamic_screen_max_length` nt — default 60, Primer3's THAL_MAX_ALIGN; opt-in up to 10 000 — the Primer3 ntthal self-dimer, 3′ self-dimer and hairpin Tm limit, default 47 °C; longer probes use Primer3's alignment-mode self_any / self_end ≤ 12.00 self-dimer screen (PRIMER_INTERNAL_MAX_SELF_ANY/_END, dpal, no length limit) and a sequence-only inverted-repeat hairpin screen) and simple repeats, and returns up to `max_probes` probes sorted by score (descending). **The additive score is a library heuristic** (its penalty values have no published source). For the **sourced ranking** set `parameters.Ranking = Primer3Penalty`: the same candidates are ranked by Primer3's internal-oligo objective `p_obj_fn` (PRIMER_INTERNAL_n_PENALTY = |Tm − `OptTm`| + |length − `OptLength`| with Primer3's default internal-oligo weights; `OptTm` / `OptLength` default to Primer3's 60 °C / 20 nt and must lie within the Tm / length window — e.g. Microarray needs its own optima — otherwise the call fails like Primer3 `_pr_data_control`) in Primer3 `primer_rec_comp` order (penalty ascending, start descending, length ascending); each probe then carries `primer3Penalty` (= primer3-py `design_primers` pick_hyb_probe_only PRIMER_INTERNAL_n_PENALTY). Primer3's complete picker is [design_probes_primer3](design_probes_primer3.md). Use a `ProbeParameters` preset (`Microarray` default, `FISH`, `NorthernBlot`, `qPCR` — Tm 68–70 °C, GC 30–80 % per the Applied Biosystems TaqMan probe guidelines —, `SouthernBlot`) or custom values. A target shorter than the minimum probe length returns an empty list.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L843](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L843)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_sequence` | string | Yes | Target DNA sequence, non-empty. |
| `parameters` | object | No | Optional `ProbeParameters` (defaults to Microarray); `Ranking` = `AdditiveScore` (default, library heuristic) or `Primer3Penalty` (Primer3 internal-oligo `p_obj_fn` around `OptTm` / `OptLength`, defaults 60 / 20). |
| `max_probes` | integer | No | Maximum probes to return (> 0, default 10). |
| `thermodynamic_screen_max_length` | integer | No | THAL_MAX_ALIGN of the ntthal self-structure screen (`ProbeParameters.ThermodynamicScreenMaxLength`): A/C/G/T probes up to this length get the ntthal screen. Default null = the parameters' value (60 = Primer3); opt-in 61–10 000 screens longer probes with the unchanged ntthal recursions (= thal.c compiled with `-DTHAL_MAX_ALIGN`; cost O(n²) per probe). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | array | Up to `max_probes` `Probe` records (sequence, start, end, Tm, GC, score, type, warnings, primer3Penalty), score-descending (default) or Primer3-penalty-ascending (`Ranking = Primer3Penalty`; `primer3Penalty` null otherwise). |

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
