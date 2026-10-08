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

Scores a single primer and returns a candidate record: length, GC%, Tm (Primer3 SantaLucia 1998 NN Tm at `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc`, Primer3 defaults 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM), longest homopolymer, the Primer3 thermodynamic secondary-structure Tm values (`hairpinTh` = primer3 `calc_hairpin` Tm, `selfAnyTh` = `calc_homodimer` Tm, `selfEndTh` = `calc_end_stability(p, p)` Tm at the same conditions; 0 when no structure or Tm < 0; the ntthal engines are bit-exact ports of primer3-py 2.3.1 `thal.c`), `hasHairpin` (= `hairpinTh` > `MaxStructureTm`, default 47 °C = Primer3 `PRIMER_MAX_HAIRPIN_TH`; with `StructureScreen = Heuristic` it is the sequence-only `HasHairpinPotential` and the `*Th` fields are null), 3′-end ΔG°37 stability, a list of quality issues (against the supplied or default `PrimerParameters`), an overall validity flag, an informational numeric score and the Primer3 per-primer `penalty` (|Tm − OptimalTm| + |length − OptimalLength|; lower is better). `position` and `is_forward` are informational and echoed back.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L1080](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L1080)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Primer sequence (non-empty). |
| `position` | integer | Yes | 0-based location (informational). |
| `is_forward` | boolean | Yes | Forward/reverse flag. |
| `parameters` | object | No | Optional design parameters. |
| `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc` | number | No | PRIMER_SALT_MONOVALENT (mM, > 0; default 50) / PRIMER_SALT_DIVALENT (mM, ≥ 0; 1.5) / PRIMER_DNTP_CONC (mM, ≥ 0; 0.6) / PRIMER_DNA_CONC (nM, > 0; 50) for the Tm and ntthal values; illegal values → `ArgumentOutOfRangeException`. |
| `opt_gc_percent` / `wt_gc_percent_gt` / `wt_gc_percent_lt` | number | No | PRIMER_OPT_GC_PERCENT (default undefined, as in Primer3's code) and PRIMER_WT_GC_PERCENT_GT / _LT (default 0) of the `penalty`; a non-zero GC weight without `opt_gc_percent` is Primer3's error "Primer GC content is part of objective function while optimum gc_content is not defined". |
| `gc_clamp` / `max_end_gc` / `max_end_stability` | integer / integer / number | No | Primer3 3′-end checks: PRIMER_GC_CLAMP (default 0), PRIMER_MAX_END_GC (0–5, default 5), PRIMER_MAX_END_STABILITY (kcal/mol ≥ 0, default 100); a failure adds an issue naming the Primer3 tag (primer3-py `check_primers` parity). |
| `annealing_temp` / `min_bound` / `max_bound` / `opt_bound` / `wt_bound_gt` / `wt_bound_lt` | number | No | PRIMER_ANNEALING_TEMP (°C ≤ 100; default −10 = off) and the fraction-bound settings PRIMER_MIN/MAX/OPT_BOUND (−10 / 110 / 97 %) and PRIMER_WT_BOUND_GT / _LT (0). With `annealing_temp` > 0 the primer's fraction bound (Primer3 `oligotm`) is reported as `bound`, a value outside [min, max] adds an issue naming PRIMER_MIN_BOUND / PRIMER_MAX_BOUND, and the weights add w × |bound − opt| to the penalty (primer3-py `check_primers` parity, F44). `opt_bound` outside [min, max] or `annealing_temp` > 100 → `ArgumentOutOfRangeException`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` / `gcContent` / `meltingTemperature` | number | Basic metrics. |
| `homopolymerLength` / `hasHairpin` / `stability3Prime` | mixed | Structural metrics. |
| `selfAnyTh` / `selfEndTh` / `hairpinTh` | number \| null | Primer3 `PRIMER_*_SELF_ANY_TH` / `_SELF_END_TH` / `_HAIRPIN_TH` (°C). |
| `selfAny` / `selfEnd` | number \| null | Primer3 alignment-mode `PRIMER_*_SELF_ANY` / `_SELF_END` (dpal scores) when `parameters.StructureScreen = Primer3Alignment` (limits `MaxSelfAny` 8 / `MaxSelfEnd` 3); null otherwise. |
| `isValid` / `issues` / `score` / `penalty` | mixed | QC verdict; `penalty` is the Primer3 ranking penalty. |
| `bound` | number \| null | Primer3 `PRIMER_LEFT/RIGHT_0_BOUND` (% bound at `annealing_temp`) when `annealing_temp` > 0; null otherwise. |
| `libraryMispriming` / `libraryMisprimingName` / `templateMispriming` / `positionPenalty` / `minSequenceQuality` / `maskFailureRate` | number / string / integer \| null | Filled only by `design_primers` (library, template, position-penalty, sequence-quality and masking options); always null here (no template). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATCGATCGATCGATCGATCG`

20-mer, 50% GC, Primer3-default Tm = 57.36 °C (primer3-py `calc_tm`) → `57.4`, penalty 2.637, homopolymer length 1.

## See Also

- [design_primers](design_primers.md), [primer_melting_temperature](primer_melting_temperature.md), [three_prime_stability](three_prime_stability.md)
