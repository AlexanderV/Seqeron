# design_probes_primer3

Pick hybridization probes exactly as Primer3 does for `PRIMER_TASK=pick_hyb_probe_only`.

## Overview

| Property | Value |
|----------|-------|
| **Server** | OligoThermodynamics |
| **Tool Name** | `design_probes_primer3` |
| **Method ID** | `ThermodynamicProbeDesigner.DesignProbesPrimer3` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of Primer3's internal-oligo (hybridization-probe) picker (`libprimer3.cc` `make_internal_oligo_list` → `pick_primer_range` → `calc_and_check_oligo_features` → `p_obj_fn`, default thermodynamic mode). Every window of `min_size..max_size` bases made only of A/C/G/T is accepted when its G+C % is within `[min_gc_percent, max_gc_percent]`, its longest mononucleotide run ≤ `max_poly_x`, its Tm (Primer3 `seqtm`, SantaLucia 1998 nearest-neighbour, at the stated salt / oligo conditions) within `[min_tm, max_tm]`, and its ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ `max_self_any_th` / `max_self_end_th` / `max_hairpin_th`. As in Primer3, windows are enumerated per 3′ end with growing length and an N, a poly-X run or a too-stable self-dimer stops the 5′ extension. Accepted probes are ranked by the Primer3 penalty `|Tm − opt_tm| + |length − opt_size|` and ordered penalty ascending, then start descending, then length ascending (`primer_rec_comp`). Verified against primer3-py 2.3.1 `design_primers` (positions, order, Tm, penalty, SELF_ANY_TH, SELF_END_TH, HAIRPIN_TH identical on 950 random runs). Defaults are Primer3's `PRIMER_INTERNAL_*` defaults.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OligoThermodynamics/ThermodynamicProbeDesigner.Core.cs#L1645](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OligoThermodynamics/ThermodynamicProbeDesigner.Core.cs#L1645)
- Algorithm doc: [Hybridization_Probe_Design.md](../../../algorithms/MolTools/Hybridization_Probe_Design.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | Template DNA sequence (probes are picked on this strand; case-insensitive). |
| `num_return` | integer | No | PRIMER_NUM_RETURN: maximum probes to return (default 5). |
| `min_size` | integer | No | PRIMER_INTERNAL_MIN_SIZE (default 18). |
| `opt_size` | integer | No | PRIMER_INTERNAL_OPT_SIZE (default 20). |
| `max_size` | integer | No | PRIMER_INTERNAL_MAX_SIZE (default 27; at most 36). |
| `min_tm` | number | No | PRIMER_INTERNAL_MIN_TM in °C (default 57). |
| `opt_tm` | number | No | PRIMER_INTERNAL_OPT_TM in °C (default 60). |
| `max_tm` | number | No | PRIMER_INTERNAL_MAX_TM in °C (default 63). |
| `min_gc_percent` | number | No | PRIMER_INTERNAL_MIN_GC in percent (default 20). |
| `max_gc_percent` | number | No | PRIMER_INTERNAL_MAX_GC in percent (default 80). |
| `max_poly_x` | integer | No | PRIMER_INTERNAL_MAX_POLY_X (default 5). |
| `max_self_any_th` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY_TH in °C (default 47). |
| `max_self_end_th` | number | No | PRIMER_INTERNAL_MAX_SELF_END_TH in °C (default 47). |
| `max_hairpin_th` | number | No | PRIMER_INTERNAL_MAX_HAIRPIN_TH in °C (default 47). |
| `monovalent_mm` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT in mM (default 50). |
| `divalent_mm` | number | No | PRIMER_INTERNAL_SALT_DIVALENT (Mg²⁺) in mM (default 0). |
| `dntp_mm` | number | No | PRIMER_INTERNAL_DNTP_CONC in mM (default 0). |
| `dna_conc_nm` | number | No | PRIMER_INTERNAL_DNA_CONC in nM (default 50). |
| `thermodynamic_oligo_alignment` | boolean | No | PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: true (default) = ntthal Tm limits; false = Primer3 alignment-mode dpal self_any / self_end limits. |
| `max_self_any` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY (alignment mode; default 12). (default: 12) |
| `max_self_end` | number | No | PRIMER_INTERNAL_MAX_SELF_END (alignment mode; default 12). (default: 12) |
| `mishyb_library` | object | No | PRIMER_INTERNAL_MISHYB_LIBRARY as a name -> sequence object (primer3-py mishyb_lib), e.g. {"Alu*2": "GGCCGGGCGCGG..."}; an optional '*weight' (0-100) after the name scales that entry; IUPAC codes allowed. Each probe is scored against every entry and its reverse complement (Primer3 dpal, unanchored local alignment); reported as libraryMishyb/libraryMishybName. |
| `max_library_mishyb` | number | No | PRIMER_INTERNAL_MAX_LIBRARY_MISHYB: maximum weighted library score of a probe (default 12). |
| `wt_library_mishyb` | number | No | PRIMER_INTERNAL_WT_LIBRARY_MISHYB: penalty weight of the library score (default 0; needs mishyb_library). |
| `lib_ambiguity_codes_consensus` | boolean | No | PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS: false (Primer3 default 0) = IUPAC codes in the library never match; true = they match every base they represent. (default: false) |
| `sequence_quality` | integer[] | No | SEQUENCE_QUALITY: one integer base quality per template base (length = template length, values within [quality_range_min, quality_range_max]); each probe then reports minSequenceQuality. |
| `min_quality` | integer | No | PRIMER_INTERNAL_MIN_QUALITY: minimum base quality of a probe (default 0; non-zero needs sequence_quality). |
| `quality_range_min` | integer | No | PRIMER_QUALITY_RANGE_MIN (default 0). |
| `quality_range_max` | integer | No | PRIMER_QUALITY_RANGE_MAX (default 100); the quality penalty is wt_seq_qual x (quality_range_max - min quality). |
| `wt_seq_qual` | number | No | PRIMER_INTERNAL_WT_SEQ_QUAL: penalty weight of (quality_range_max - the probe's minimum base quality) (default 0; needs sequence_quality). |
| `wt_end_qual` | number | No | PRIMER_INTERNAL_WT_END_QUAL: accepted for Primer3 compatibility; no effect (as in Primer3 2.3.1). (default: 0) |
| `opt_gc_percent` | number | No | PRIMER_INTERNAL_OPT_GC_PERCENT: GC optimum of the GC penalty terms (default undefined, as in Primer3's code; required when wt_gc_percent_gt/lt != 0 - Primer3 'Hyb probe GC content is part of objective function while optimum gc_content is not defined'). |
| `wt_gc_percent_gt` | number | No | PRIMER_INTERNAL_WT_GC_PERCENT_GT: penalty weight for GC% above opt_gc_percent (default 0). |
| `wt_gc_percent_lt` | number | No | PRIMER_INTERNAL_WT_GC_PERCENT_LT: penalty weight for GC% below opt_gc_percent (default 0). |
| `lowercase_masking` | boolean | No | PRIMER_LOWERCASE_MASKING: reject probes whose 3'-terminal template base is lower case (a/c/g/t of the template as given; lower case elsewhere is accepted). Default false. (default: false) |
| `annealing_temp` | number | No | PRIMER_ANNEALING_TEMP: annealing temperature in °C (<= 100; default -10 = off). When > 0 each probe's fraction bound at this temperature (Primer3 oligotm, at the probe conditions) is reported as bound and probes outside [min_bound, max_bound] are rejected. (default: -10) |
| `min_bound` | number | No | PRIMER_INTERNAL_MIN_BOUND: minimum fraction bound, % (default -10; only with annealing_temp > 0). |
| `max_bound` | number | No | PRIMER_INTERNAL_MAX_BOUND: maximum fraction bound, % (default 110; only with annealing_temp > 0). |
| `opt_bound` | number | No | PRIMER_INTERNAL_OPT_BOUND: optimum of the bound penalty terms, % (default 97; must lie in [min_bound, max_bound]). |
| `wt_bound_gt` | number | No | PRIMER_INTERNAL_WT_BOUND_GT: penalty weight of the fraction bound above opt_bound (default 0). |
| `wt_bound_lt` | number | No | PRIMER_INTERNAL_WT_BOUND_LT: penalty weight of the fraction bound below opt_bound (default 0). As in Primer3 the internal-oligo bound terms are not gated by annealing_temp: without it the bound is OLIGOTM_ERROR -999999.9999, so wt_bound_lt adds wt_bound_lt x (opt_bound + 999999.9999). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | object[] | Best first. Each: `sequence`, `start` (0-based), `length`, `tm`, `gcPercent`, `selfAnyTh`, `selfEndTh`, `hairpinTh`, `penalty` (the Primer3 `PRIMER_INTERNAL_n_*` values); with `mishyb_library`: `libraryMishyb` / `libraryMishybName` (PRIMER_INTERNAL_n_LIBRARY_MISHYB score and entry; primer3-py key `…_LIBRARY_MISPRIMING`); with `sequence_quality`: `minSequenceQuality`; with `annealing_temp` > 0: `bound` (PRIMER_INTERNAL_n_BOUND, %). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Template sequence cannot be null or empty |
| 1002 | num_return must be at least 1 (Primer3: PRIMER_NUM_RETURN < 1) |
| 1004 | Internal oligo mispriming score is part of objective function while mishyb library is not defined |
| 1003 | Sizes must satisfy 1 ≤ min_size ≤ max_size ≤ 36 |

## Examples

### Example 1: Primer3 defaults

120-nt template `CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC`.

**Response (abridged):** first probe `start = 27`, `length = 22`, `tm = 57.5737`, `penalty = 4.4263`, `selfAnyTh = 15.4317`, `selfEndTh = 1.5232`, `hairpinTh = 37.4725` (= primer3-py `PRIMER_INTERNAL_0_*`); positions of the 5 probes (27,22), (27,23), (67,24), (68,23), (67,23).

## See Also

- [design_probes](design_probes.md), [validate_probe](validate_probe.md), [design_primers](design_primers.md)
