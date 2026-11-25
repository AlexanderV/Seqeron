# evaluate_primer

Evaluate a single primer against quality criteria.

## Overview

| Property | Value |
|----------|-------|
| **Server** | OligoThermodynamics |
| **Tool Name** | `evaluate_primer` |
| **Method ID** | `ThermodynamicPrimerDesigner.EvaluatePrimer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scores a single primer and returns a candidate record: length, GC%, Tm (Primer3 SantaLucia 1998 NN Tm at `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc`, Primer3 defaults 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM), longest homopolymer, the Primer3 thermodynamic secondary-structure Tm values (`hairpinTh` = primer3 `calc_hairpin` Tm, `selfAnyTh` = `calc_homodimer` Tm, `selfEndTh` = `calc_end_stability(p, p)` Tm at the same conditions; 0 when no structure or Tm < 0; the ntthal engines are bit-exact ports of primer3-py 2.3.1 `thal.c`), `hasHairpin` (= `hairpinTh` > `MaxStructureTm`, default 47 °C = Primer3 `PRIMER_MAX_HAIRPIN_TH`; with `StructureScreen = Heuristic` it is the sequence-only `HasHairpinPotential` and the `*Th` fields are null), 3′-end ΔG°37 stability, a list of quality issues (against the supplied or default `PrimerParameters`), an overall validity flag, an informational numeric score and the Primer3 per-primer `penalty` (|Tm − OptimalTm| + |length − OptimalLength|; lower is better). `position` and `is_forward` are informational and echoed back.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L1876](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L1876)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Primer sequence to evaluate. |
| `position` | integer | Yes | 0-based location of the primer in the template (informational). |
| `is_forward` | boolean | Yes | True if this is a forward primer; false for reverse. |
| `parameters` | object | No | Optional primer design parameters. |
| `salt_monovalent` | number | No | PRIMER_SALT_MONOVALENT: monovalent cation concentration in mM (> 0) for the Tm and ntthal structure values (default 50). |
| `salt_divalent` | number | No | PRIMER_SALT_DIVALENT: Mg2+ concentration in mM (>= 0; default 1.5). |
| `dntp_conc` | number | No | PRIMER_DNTP_CONC: dNTP concentration in mM (>= 0; default 0.6). |
| `dna_conc` | number | No | PRIMER_DNA_CONC: primer concentration in nM (> 0; default 50). |
| `opt_gc_percent` | number | No | PRIMER_OPT_GC_PERCENT: GC optimum of the GC penalty terms (default undefined, as in Primer3's code; required when wt_gc_percent_gt/lt != 0 - Primer3 'Primer GC content is part of objective function while optimum gc_content is not defined'). |
| `wt_gc_percent_gt` | number | No | PRIMER_WT_GC_PERCENT_GT: penalty weight for GC% above opt_gc_percent (default 0). |
| `wt_gc_percent_lt` | number | No | PRIMER_WT_GC_PERCENT_LT: penalty weight for GC% below opt_gc_percent (default 0). |
| `gc_clamp` | integer | No | PRIMER_GC_CLAMP: number of consecutive G/C required at the 3' end (default 0). |
| `max_end_gc` | integer | No | PRIMER_MAX_END_GC: maximum G/C among the five 3'-most bases (0-5; default 5 = no limit). |
| `max_end_stability` | number | No | PRIMER_MAX_END_STABILITY: maximum 3' end stability (-dG of the 3' pentamer, kcal/mol, >= 0; default 100 = no limit). |
| `annealing_temp` | number | No | PRIMER_ANNEALING_TEMP: annealing temperature in °C (<= 100; default -10 = off). When > 0 the fraction bound at this temperature is reported as bound and checked against min_bound / max_bound. |
| `min_bound` | number | No | PRIMER_MIN_BOUND: minimum fraction bound, % (default -10; only with annealing_temp > 0). |
| `max_bound` | number | No | PRIMER_MAX_BOUND: maximum fraction bound, % (default 110; only with annealing_temp > 0). |
| `opt_bound` | number | No | PRIMER_OPT_BOUND: optimum fraction bound of the bound penalty terms, % (default 97). |
| `wt_bound_gt` | number | No | PRIMER_WT_BOUND_GT: penalty weight of the fraction bound above opt_bound (default 0; only with annealing_temp > 0). |
| `wt_bound_lt` | number | No | PRIMER_WT_BOUND_LT: penalty weight of the fraction bound below opt_bound (default 0; only with annealing_temp > 0). |

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
