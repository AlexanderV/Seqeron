# design_probes_primer3

Pick hybridization probes exactly as Primer3 does for `PRIMER_TASK=pick_hyb_probe_only`.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_probes_primer3` |
| **Method ID** | `ProbeDesigner.DesignProbesPrimer3` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of Primer3's internal-oligo (hybridization-probe) picker (`libprimer3.cc` `make_internal_oligo_list` → `pick_primer_range` → `calc_and_check_oligo_features` → `p_obj_fn`, default thermodynamic mode). Every window of `min_size..max_size` bases made only of A/C/G/T is accepted when its G+C % is within `[min_gc_percent, max_gc_percent]`, its longest mononucleotide run ≤ `max_poly_x`, its Tm (Primer3 `seqtm`, SantaLucia 1998 nearest-neighbour, at the stated salt / oligo conditions) within `[min_tm, max_tm]`, and its ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ `max_self_any_th` / `max_self_end_th` / `max_hairpin_th`. As in Primer3, windows are enumerated per 3′ end with growing length and an N, a poly-X run or a too-stable self-dimer stops the 5′ extension. Accepted probes are ranked by the Primer3 penalty `|Tm − opt_tm| + |length − opt_size|` and ordered penalty ascending, then start descending, then length ascending (`primer_rec_comp`). Verified against primer3-py 2.3.1 `design_primers` (positions, order, Tm, penalty, SELF_ANY_TH, SELF_END_TH, HAIRPIN_TH identical on 950 random runs). Defaults are Primer3's `PRIMER_INTERNAL_*` defaults.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1672](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1672)
- Algorithm doc: [Hybridization_Probe_Design.md](../../../algorithms/MolTools/Hybridization_Probe_Design.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | Template DNA (non-empty; case-insensitive unless `lowercase_masking`; probes are picked on this strand). |
| `num_return` | integer | No | PRIMER_NUM_RETURN (≥ 1, default 5; Primer3 rejects < 1). |
| `min_size` / `opt_size` / `max_size` | integer | No | PRIMER_INTERNAL_MIN/OPT/MAX_SIZE (18 / 20 / 27; 1 ≤ min ≤ max ≤ 36). |
| `min_tm` / `opt_tm` / `max_tm` | number | No | PRIMER_INTERNAL_MIN/OPT/MAX_TM in °C (57 / 60 / 63). |
| `min_gc_percent` / `max_gc_percent` | number | No | PRIMER_INTERNAL_MIN/MAX_GC (20 / 80 %). |
| `max_poly_x` | integer | No | PRIMER_INTERNAL_MAX_POLY_X (5). |
| `max_self_any_th` / `max_self_end_th` / `max_hairpin_th` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH (47 °C). |
| `monovalent_mm` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT (50 mM). |
| `divalent_mm` | number | No | PRIMER_INTERNAL_SALT_DIVALENT (0 mM). |
| `dntp_mm` | number | No | PRIMER_INTERNAL_DNTP_CONC (0 mM). |
| `dna_conc_nm` | number | No | PRIMER_INTERNAL_DNA_CONC (50 nM). |
| `thermodynamic_oligo_alignment` | boolean | No | PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT (default true). false = Primer3 alignment mode: dpal `self_any` ≤ `max_self_any`, `self_end` ≤ `max_self_end`, no hairpin; `selfAny`/`selfEnd` reported, `*Th` fields NaN. |
| `max_self_any` / `max_self_end` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END (12, alignment mode). |
| `mishyb_library` | object | No | PRIMER_INTERNAL_MISHYB_LIBRARY as a name → sequence object (primer3-py `mishyb_lib`), e.g. `{"Alu*2": "GGCCGGGCGCGG…"}`; same entry format as `design_primers` `mispriming_library` (`*weight` 0–100, IUPAC codes, "reverse <name>" entries added). Each probe is aligned with every entry by Primer3's dpal unanchored local alignment (+1/−1/−0.25 N/−2 gap, max gap 1); a probe whose weighted score exceeds `max_library_mishyb` is rejected (and ends the 5′ extension of its 3′ end, as in Primer3). |
| `max_library_mishyb` | number | No | PRIMER_INTERNAL_MAX_LIBRARY_MISHYB (default 12; compared as a C `short`; > 32767 rejected in alignment mode). |
| `wt_library_mishyb` | number | No | PRIMER_INTERNAL_WT_LIBRARY_MISHYB (default 0): weight of the library score in the penalty; non-zero without `mishyb_library` → error (Primer3 `_pr_data_control`). |
| `lib_ambiguity_codes_consensus` | boolean | No | PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS (default false = Primer3 0: IUPAC codes in the library never match; true: they match every base they represent). Verified against primer3-py 2.3.1 `design_primers(mishyb_lib=…)` (B07 F42). |
| `sequence_quality` | integer[] | No | SEQUENCE_QUALITY: one integer base quality per template base (length = template length, values within [`quality_range_min`, `quality_range_max`]); each probe reports `minSequenceQuality` (PRIMER_INTERNAL_n_MIN_SEQ_QUALITY). |
| `min_quality` | integer | No | PRIMER_INTERNAL_MIN_QUALITY (default 0): probes whose minimum base quality is lower are rejected (a five-prime problem, as in Primer3); non-zero without `sequence_quality` → error ("Sequence quality data missing"). |
| `quality_range_min` / `quality_range_max` | integer | No | PRIMER_QUALITY_RANGE_MIN / _MAX (defaults 0 / 100). |
| `wt_seq_qual` / `wt_end_qual` | number | No | PRIMER_INTERNAL_WT_SEQ_QUAL (default 0): weight × (`quality_range_max` − minimum base quality), needs `sequence_quality`; PRIMER_INTERNAL_WT_END_QUAL is accepted with no effect (never read by Primer3 2.3.1). Verified against primer3-py 2.3.1 (B07 F45). |
| `opt_gc_percent` | number | No | PRIMER_INTERNAL_OPT_GC_PERCENT (default undefined, as in Primer3's code); required when `wt_gc_percent_gt` / `wt_gc_percent_lt` ≠ 0 (Primer3 "Hyb probe GC content is part of objective function while optimum gc_content is not defined"). |
| `wt_gc_percent_gt` / `wt_gc_percent_lt` | number | No | PRIMER_INTERNAL_WT_GC_PERCENT_GT / _LT (default 0): weight × the GC% deviation above / below `opt_gc_percent`. |
| `lowercase_masking` | boolean | No | PRIMER_LOWERCASE_MASKING (default false): reject probes whose 3′-terminal template base is lower case (a/c/g/t of the template as given); lower case elsewhere is accepted. |
| `annealing_temp` | number | No | PRIMER_ANNEALING_TEMP in °C (≤ 100; default −10 = off). When > 0 each probe's fraction bound at this temperature (Primer3 `oligotm` at the probe conditions) is reported as `bound` and probes outside [`min_bound`, `max_bound`] are rejected. |
| `min_bound` / `max_bound` / `opt_bound` | number | No | PRIMER_INTERNAL_MIN_BOUND / _MAX_BOUND / _OPT_BOUND in % (−10 / 110 / 97); `opt_bound` must lie in [`min_bound`, `max_bound`] (Primer3 `_pr_data_control`). |
| `wt_bound_gt` / `wt_bound_lt` | number | No | PRIMER_INTERNAL_WT_BOUND_GT / _LT (default 0): weight × the bound deviation above / below `opt_bound`. Not gated by `annealing_temp` (as in Primer3): without it the bound is OLIGOTM_ERROR −999999.9999, so `wt_bound_lt` adds `wt_bound_lt` × (`opt_bound` + 999999.9999). Verified against primer3-py 2.3.1 (B07 F47). |

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
