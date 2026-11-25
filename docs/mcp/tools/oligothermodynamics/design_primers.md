# design_primers

Design a forward/reverse PCR primer pair flanking a target region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | OligoThermodynamics |
| **Tool Name** | `design_primers` |
| **Method ID** | `ThermodynamicPrimerDesigner.DesignPrimers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Implements Primer3's pair search (`libprimer3.cc` `choose_pair_or_triple` / `characterize_pair` / `obj_fn`; verified against primer3-py 2.3.1 `design_primers`):

- **Candidates.** Forward candidates end at or before `target_start`, reverse candidates (evaluated on the reverse complement) start at or after `target_end`; only positions that can give a product in `product_size_range` are enumerated. A candidate is kept when it passes every per-primer limit of `PrimerParameters` (length, GC %, Primer3 SantaLucia 1998 NN Tm at `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc` — Primer3 defaults 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM; the same conditions drive the ntthal structure / pair complementarity values and the product Tm, as in Primer3 — poly-X, dinucleotide repeat, no non-ACGT base).
- **Real filters of a pair** (in Primer3's order): product size in the current `product_size_range` range (default **100–300 bp**, Primer3 PRIMER_PRODUCT_SIZE_RANGE; ranges are tried in order), `|Tm_f − Tm_r| ≤ max_tm_difference` (default **5 °C**; Primer3's own default is 100), and — with the default `PrimerStructureScreen.Primer3Thermodynamic` — Primer3's thermodynamic structure screen: each primer's ntthal self-dimer (SELF_ANY_TH), 3′ self-dimer (SELF_END_TH) and hairpin (HAIRPIN_TH) Tm ≤ **47 °C**, and the pair's ntthal hetero-dimer (COMPL_ANY_TH) and 3′ hetero-dimer (COMPL_END_TH) Tm ≤ **47 °C** (`MaxStructureTm`). With `StructureScreen = Heuristic` the screen is the sequence-only stem-loop test plus the alignment 3′ complementarity test. With `pick_internal_oligo` the pair also needs an internal hybridization oligo strictly between the primers (Primer3 PRIMER_INTERNAL_* defaults).
- **Selection.** The pair with the lowest Primer3 pair penalty (PRIMER_PAIR_WT_PR_PENALTY = 1: sum of the per-primer penalties |Tm − OptimalTm| + |length − OptimalLength|); ties as Primer3's `compare_primer_pair`. `pairs` lists up to `num_return` ranked pairs (PRIMER_NUM_RETURN; each later pair excludes the earlier ones). When no pair qualifies, the individually best primers are returned with `isValid = false` and a message naming the violated constraint. `product_size = reverse.Position + reverse.Length − forward.Position`.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L2113](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L2113)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | DNA template (A/C/G/T). |
| `target_start` | integer | Yes | 0-based inclusive start of target region. |
| `target_end` | integer | Yes | 0-based exclusive end of target region (with the default inside_penalty/outside_penalty primers never overlap [target_start, target_end)). |
| `parameters` | object | No | Optional primer design parameters (lengths, GC%, Tm, repeats, GC-clamp/3' stability checks, structure screen). Defaults are used if null. |
| `product_size_range` | string | No | PRIMER_PRODUCT_SIZE_RANGE in Primer3 syntax, e.g. "100-300" or "150-250 100-400" (ranges in order of preference; default 100-300). |
| `max_tm_difference` | number | No | PRIMER_PAIR_MAX_DIFF_TM: maximum \|Tm_forward - Tm_reverse\| in °C (default 5; Primer3's own default is 100). |
| `num_return` | integer | No | PRIMER_NUM_RETURN: number of ranked pairs listed in 'pairs' (default 1). |
| `pick_internal_oligo` | boolean | No | PRIMER_PICK_INTERNAL_OLIGO: also pick an internal hybridization oligo between the primers (Primer3 PRIMER_INTERNAL_* defaults; default false). (default: false) |
| `pair_max_compl_any` | number | No | PRIMER_PAIR_MAX_COMPL_ANY: maximum Primer3 alignment-mode pair compl_any (used only with parameters.StructureScreen = Primer3Alignment; default 8). (default: 8) |
| `pair_max_compl_end` | number | No | PRIMER_PAIR_MAX_COMPL_END: maximum Primer3 alignment-mode pair compl_end (used only with parameters.StructureScreen = Primer3Alignment; default 3). (default: 3) |
| `salt_monovalent` | number | No | PRIMER_SALT_MONOVALENT: monovalent cation concentration in mM (> 0) for primer Tm, ntthal structure / pair complementarity and product Tm (default 50). |
| `salt_divalent` | number | No | PRIMER_SALT_DIVALENT: Mg2+ concentration in mM (>= 0; default 1.5). |
| `dntp_conc` | number | No | PRIMER_DNTP_CONC: dNTP concentration in mM (>= 0; default 0.6). |
| `dna_conc` | number | No | PRIMER_DNA_CONC: primer concentration in nM (> 0; default 50). |
| `opt_gc_percent` | number | No | PRIMER_OPT_GC_PERCENT: GC optimum of the primer GC penalty terms (default undefined, as in Primer3's code; required when wt_gc_percent_gt/lt != 0 - otherwise Primer3's error 'Primer GC content is part of objective function while optimum gc_content is not defined'). |
| `wt_gc_percent_gt` | number | No | PRIMER_WT_GC_PERCENT_GT: penalty weight for GC% above opt_gc_percent (default 0). |
| `wt_gc_percent_lt` | number | No | PRIMER_WT_GC_PERCENT_LT: penalty weight for GC% below opt_gc_percent (default 0). |
| `internal_salt_monovalent` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT for the internal oligo, mM (> 0; default 50). |
| `internal_salt_divalent` | number | No | PRIMER_INTERNAL_SALT_DIVALENT for the internal oligo, mM (>= 0; default 0). |
| `internal_dntp_conc` | number | No | PRIMER_INTERNAL_DNTP_CONC for the internal oligo, mM (>= 0; default 0). |
| `internal_dna_conc` | number | No | PRIMER_INTERNAL_DNA_CONC for the internal oligo, nM (> 0; default 50). |
| `internal_opt_gc_percent` | number | No | PRIMER_INTERNAL_OPT_GC_PERCENT: GC optimum of the internal-oligo GC penalty terms (default undefined; required when internal_wt_gc_percent_gt/lt != 0, even without pick_internal_oligo - Primer3 'Hyb probe GC content is part of objective function while optimum gc_content is not defined'). |
| `internal_wt_gc_percent_gt` | number | No | PRIMER_INTERNAL_WT_GC_PERCENT_GT (default 0). |
| `internal_wt_gc_percent_lt` | number | No | PRIMER_INTERNAL_WT_GC_PERCENT_LT (default 0). |
| `gc_clamp` | integer | No | PRIMER_GC_CLAMP: number of consecutive G/C required at each primer's 3' end (default 0; must be <= the minimum primer length). |
| `max_end_gc` | integer | No | PRIMER_MAX_END_GC: maximum G/C among the five 3'-most bases of each primer (0-5; default 5 = no limit). |
| `max_end_stability` | number | No | PRIMER_MAX_END_STABILITY: maximum 3' end stability (Primer3 end_stability = -dG of the 3' pentamer, kcal/mol, >= 0; default 100 = no limit). |
| `min_left_three_prime_distance` | integer | No | PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: after a pair is selected, later pairs may not use a left primer whose 3' end is closer than this (0 = not the identical primer; default -1 = reuse allowed). |
| `min_right_three_prime_distance` | integer | No | PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE: as min_left_three_prime_distance for right primers (default -1). |
| `min_three_prime_distance` | integer | No | PRIMER_MIN_THREE_PRIME_DISTANCE: sets both the left and right minimum 3' distances (cannot be combined with min_left/right_three_prime_distance). |
| `mispriming_library` | object | No | PRIMER_MISPRIMING_LIBRARY as a name -> sequence object (primer3-py misprime_lib), e.g. {"Alu*2": "GGCCGGGCGCGG..."}; an optional '*weight' (0-100) after the name scales that entry; IUPAC codes allowed. Each primer is scored against every entry and its reverse complement (Primer3 dpal, 3'-anchored); reported as libraryMispriming/libraryMisprimingName per primer and pair. |
| `max_library_mispriming` | number | No | PRIMER_MAX_LIBRARY_MISPRIMING: maximum weighted library score of one primer (default 12). |
| `pair_max_library_mispriming` | number | No | PRIMER_PAIR_MAX_LIBRARY_MISPRIMING: maximum pair library score (integer part of left + right score for the same entry; default 24). |
| `wt_library_mispriming` | number | No | PRIMER_WT_LIBRARY_MISPRIMING: per-primer penalty weight of the library score (default 0; needs mispriming_library). |
| `pair_wt_library_mispriming` | number | No | PRIMER_PAIR_WT_LIBRARY_MISPRIMING: pair penalty weight of the pair library score (default 0; needs mispriming_library). |
| `lib_ambiguity_codes_consensus` | boolean | No | PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS: false (Primer3 default 0) = IUPAC codes in the library never match; true = they match every base they represent (N matches anything). |
| `internal_mishyb_library` | object | No | PRIMER_INTERNAL_MISHYB_LIBRARY as a name -> sequence object (primer3-py mishyb_lib; same format as mispriming_library), used with pick_internal_oligo: each internal oligo is scored against every entry and its reverse complement (Primer3 dpal, unanchored local alignment); reported as internalOligo.libraryMishyb/libraryMishybName. |
| `internal_max_library_mishyb` | number | No | PRIMER_INTERNAL_MAX_LIBRARY_MISHYB: maximum weighted library score of the internal oligo (default 12). |
| `internal_wt_library_mishyb` | number | No | PRIMER_INTERNAL_WT_LIBRARY_MISHYB: internal-oligo penalty weight of the library score (default 0; needs internal_mishyb_library). |
| `thermodynamic_template_alignment` | boolean | No | PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT: false (Primer3 default 0) = template mispriming is a dpal score (3'-anchored local alignment of the primer with the template outside its own site and with the opposite strand); true = the ntthal 3'-end (END1) Tm in °C (template <= 10000 nt). Uses the *_template_mispriming_th limits/weights. |
| `max_template_mispriming` | number | No | PRIMER_MAX_TEMPLATE_MISPRIMING: maximum template mispriming score of one primer (alignment mode; default -100 = not checked); reported as templateMispriming per primer. |
| `max_template_mispriming_th` | number | No | PRIMER_MAX_TEMPLATE_MISPRIMING_TH: maximum template mispriming Tm of one primer, °C (thermodynamic mode; default -100 = not checked). |
| `pair_max_template_mispriming` | number | No | PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING: maximum pair template mispriming score (max(left same-strand + right other-strand, left other-strand + right same-strand); alignment mode; default -100 = not checked). |
| `pair_max_template_mispriming_th` | number | No | PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH: maximum pair template mispriming value (thermodynamic mode; default -100). As in Primer3 a pair fails when this is non-zero and exceeded, so with the default and a non-zero pair_wt_template_mispriming_th every pair fails; 0 = no limit. |
| `wt_template_mispriming` | number | No | PRIMER_WT_TEMPLATE_MISPRIMING: per-primer penalty weight of the template mispriming score (alignment mode; default 0). |
| `wt_template_mispriming_th` | number | No | PRIMER_WT_TEMPLATE_MISPRIMING_TH: per-primer penalty weight of the template mispriming Tm (thermodynamic mode, Primer3 temp_cutoff rule; default 0). |
| `pair_wt_template_mispriming` | number | No | PRIMER_PAIR_WT_TEMPLATE_MISPRIMING: pair penalty weight of the pair template mispriming score (alignment mode; default 0). |
| `pair_wt_template_mispriming_th` | number | No | PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH: pair penalty weight of the pair template mispriming value (thermodynamic mode; default 0). |
| `annealing_temp` | number | No | PRIMER_ANNEALING_TEMP: annealing temperature in °C (<= 100; default -10 = off). When > 0 the fraction of each primer / internal oligo bound at this temperature is computed (Primer3 oligotm, SantaLucia 1998) and the bound limits / weights below apply; it is also the internal oligo's annealing temperature. |
| `min_bound` | number | No | PRIMER_MIN_BOUND: minimum fraction bound of a primer, % (default -10; only with annealing_temp > 0). |
| `max_bound` | number | No | PRIMER_MAX_BOUND: maximum fraction bound of a primer, % (default 110; only with annealing_temp > 0). |
| `opt_bound` | number | No | PRIMER_OPT_BOUND: optimum fraction bound of the bound penalty terms, % (default 97; must lie in [min_bound, max_bound]). |
| `wt_bound_gt` | number | No | PRIMER_WT_BOUND_GT: per-primer penalty weight of the fraction bound above opt_bound (default 0; only with annealing_temp > 0). |
| `wt_bound_lt` | number | No | PRIMER_WT_BOUND_LT: per-primer penalty weight of the fraction bound below opt_bound (default 0; only with annealing_temp > 0). |
| `internal_min_bound` | number | No | PRIMER_INTERNAL_MIN_BOUND: minimum fraction bound of the internal oligo, % (default -10; only with annealing_temp > 0). |
| `internal_max_bound` | number | No | PRIMER_INTERNAL_MAX_BOUND: maximum fraction bound of the internal oligo, % (default 110; only with annealing_temp > 0). |
| `internal_opt_bound` | number | No | PRIMER_INTERNAL_OPT_BOUND: optimum fraction bound of the internal-oligo bound terms, % (default 97). |
| `internal_wt_bound_gt` | number | No | PRIMER_INTERNAL_WT_BOUND_GT: internal-oligo penalty weight of the fraction bound above internal_opt_bound (default 0). |
| `internal_wt_bound_lt` | number | No | PRIMER_INTERNAL_WT_BOUND_LT: internal-oligo penalty weight of the fraction bound below internal_opt_bound (default 0). As in Primer3 it applies even without annealing_temp, where the bound is Primer3's error value -999999.9999. |
| `inside_penalty` | number | No | PRIMER_INSIDE_PENALTY: penalty per base of a primer 3' end inside the target (default -1). Any value other than the defaults -1 / 0 of inside_penalty / outside_penalty lets primers overlap the target as long as their 3' end does not pass its far end; as in Primer3 the default -1 then makes inside positions negative, and a negative pair penalty is an error. |
| `outside_penalty` | number | No | PRIMER_OUTSIDE_PENALTY: penalty per base of distance between a primer 3' end and the target (default 0; see inside_penalty). |
| `wt_pos_penalty` | number | No | PRIMER_WT_POS_PENALTY: weight of the position penalty in the primer penalty (default 1). |
| `sequence_quality` | integer[] | No | SEQUENCE_QUALITY: one integer base quality per template base (e.g. Phred scores; length = template length, values within [quality_range_min, quality_range_max]). Each primer / internal oligo then reports minSequenceQuality (its minimum base quality). |
| `min_quality` | integer | No | PRIMER_MIN_QUALITY: minimum base quality of a primer (default 0; non-zero needs sequence_quality and must lie in the quality range). |
| `min_end_quality` | integer | No | PRIMER_MIN_END_QUALITY: minimum base quality of a primer's five 3'-most bases (default 0; only with sequence_quality). |
| `quality_range_min` | integer | No | PRIMER_QUALITY_RANGE_MIN: smallest allowed sequence_quality value (default 0). |
| `quality_range_max` | integer | No | PRIMER_QUALITY_RANGE_MAX: largest allowed sequence_quality value (default 100); the quality penalty is wt_seq_qual x (quality_range_max - min quality). |
| `wt_seq_qual` | number | No | PRIMER_WT_SEQ_QUAL: per-primer penalty weight of (quality_range_max - the primer's minimum base quality) (default 0; needs sequence_quality). |
| `wt_end_qual` | number | No | PRIMER_WT_END_QUAL: accepted for Primer3 compatibility; Primer3 2.3.1 never uses it in the penalty (no effect). |
| `internal_min_quality` | integer | No | PRIMER_INTERNAL_MIN_QUALITY: minimum base quality of the internal oligo (default 0; non-zero needs sequence_quality). |
| `internal_wt_seq_qual` | number | No | PRIMER_INTERNAL_WT_SEQ_QUAL: internal-oligo penalty weight of (quality_range_max - its minimum base quality) (default 0; needs sequence_quality). |
| `internal_wt_end_qual` | number | No | PRIMER_INTERNAL_WT_END_QUAL: accepted for Primer3 compatibility; no effect (as in Primer3 2.3.1). |
| `pair_wt_io_penalty` | number | No | PRIMER_PAIR_WT_IO_PENALTY: pair penalty weight of the internal-oligo penalty (default 0; non-zero requires pick_internal_oligo, as in Primer3). |
| `mask_template` | boolean | No | PRIMER_MASK_TEMPLATE: mask the template with Primer3's k-mer masker (needs mask_kmers_11 and mask_kmers_16): primers whose 3' base is masked on their strand are rejected and each primer gets its predicted failure rate (reported as maskFailureRate; weight wt_mask_failure_rate). Default false. |
| `mask_kmers_11` | object | No | The masker's 11-mer genome counts as a k-mer -> count object (the content of PRIMER_MASK_KMERLIST_PATH/<prefix>_11.list, e.g. a GenomeTester4 glistmaker list); k-mers absent here and as reverse complement count as 1. |
| `mask_kmers_16` | object | No | The masker's 16-mer genome counts (content of <prefix>_16.list), as mask_kmers_11. |
| `mask_failure_rate` | number | No | PRIMER_MASK_FAILURE_RATE: template positions where a primer would end with a predicted failure rate above this are masked (default 0.1; 0 = none). |
| `mask_5p_direction` | integer | No | PRIMER_MASK_5P_DIRECTION: bases masked from such a 3' end towards the primer 5' end, including it (default 1). |
| `mask_3p_direction` | integer | No | PRIMER_MASK_3P_DIRECTION: bases masked beyond such a 3' end (default 0). |
| `wt_mask_failure_rate` | number | No | PRIMER_WT_MASK_FAILURE_RATE: per-primer penalty weight of the predicted failure rate (default 0; effective only with mask_template). |
| `lowercase_masking` | boolean | No | PRIMER_LOWERCASE_MASKING: reject primers and internal oligos whose 3'-terminal template base is lower case (a/c/g/t of the template as given; lower case elsewhere in an oligo is accepted). Default false; mask_template implies it. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `forward` | object \| null | Best forward `PrimerCandidate` (PRIMER_LEFT_0; null if none valid). |
| `reverse` | object \| null | Best reverse `PrimerCandidate` (PRIMER_RIGHT_0; null if none valid). |
| `isValid` | boolean | True when a compatible pair was found. |
| `message` | string | Human-readable status. |
| `productSize` | integer | Amplicon size in bp. |
| `pairPenalty` | number \| null | PRIMER_PAIR_0_PENALTY. |
| `productTm` | number \| null | PRIMER_PAIR_0_PRODUCT_TM (Primer3 `long_seq_tm`, °C). |
| `complAnyTh` / `complEndTh` | number \| null | PRIMER_PAIR_0_COMPL_ANY_TH / _COMPL_END_TH (°C); null under the alignment screen. |
| `complAny` / `complEnd` | number \| null | PRIMER_PAIR_0_COMPL_ANY / _COMPL_END (alignment screen only; verified vs primer3-py `design_primers` with PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0). |
| `templateMispriming` | number \| null | PRIMER_PAIR_0_TEMPLATE_MISPRIMING (alignment mode) or _TH (thermodynamic mode); null unless a pair template limit / weight is set. Each `PrimerCandidate` carries PRIMER_LEFT/RIGHT_0_TEMPLATE_MISPRIMING[_TH] (`templateMispriming`) when a template-mispriming setting is active. |
| `libraryMispriming` / `libraryMisprimingName` | number / string \| null | PRIMER_PAIR_0_LIBRARY_MISPRIMING (score, entry); null without `mispriming_library`. Each `PrimerCandidate` carries its PRIMER_LEFT/RIGHT_0_LIBRARY_MISPRIMING the same way (`libraryMispriming`, `libraryMisprimingName`). |
| `bound` / `positionPenalty` (per `PrimerCandidate`) | number \| null | PRIMER_LEFT/RIGHT_0_BOUND (with `annealing_temp` > 0) and PRIMER_LEFT/RIGHT_0_POSITION_PENALTY (with non-default `inside_penalty` / `outside_penalty`); null otherwise. |
| `minSequenceQuality` (per `PrimerCandidate` and internal oligo) | integer \| null | PRIMER_LEFT/RIGHT/INTERNAL_0_MIN_SEQ_QUALITY (with `sequence_quality`); null otherwise. |
| `maskFailureRate` (per `PrimerCandidate`) | number \| null | The primer's predicted failure rate (Primer3 `failure_rate`) with `mask_template`; null otherwise. |
| `internalOligo` | object \| null | PRIMER_INTERNAL_0_* (`sequence`, `start`, `length`, `tm`, `gcPercent`, `selfAnyTh`, `selfEndTh`, `hairpinTh`, `penalty`, `bound` with `annealing_temp` > 0). |
| `pairs` | array | Ranked valid pairs (k = 0 … num_return − 1), each with the fields above; empty when none qualifies. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Template cannot be null or empty |
| 1002 | Target start must be non-negative |
| 1003 | Target end must be within the template |
| 1004 | Target start must be strictly less than target end |
| 1005 | num_return must be at least 1 |
| 1006 | max_tm_difference must be non-negative |
| 1007 | Invalid product size range |
| 1008 | Optimum primer / internal oligo fraction binding lower than minimum or higher than maximum; Annealing temperature higher than 100 C (Primer3 `_pr_data_control`, `ArgumentOutOfRangeException`) |
| 1009 | Negative primer-pair penalty (Primer3 `obj_fn` `PR_ASSERT(sum >= 0.0)`, `InvalidOperationException`) |

## Examples

### Example 1: Standard 258 bp template, target [100, 150)

Forward region of `GAACTCGT` units, a 50 bp poly-T target, and a reverse region of `TCCGAAGT` units. primer3-py returns the same pair (PRIMER_LEFT_0 = [76,20], PRIMER_RIGHT_0 = [181,20], PRIMER_PAIR_0_PENALTY = 0.6983187292252637, PRIMER_PAIR_0_PRODUCT_SIZE = 106, PRIMER_PAIR_0_PRODUCT_TM = 74.016).

**Input:** `{ "target_start": 100, "target_end": 150 }`

**Response (abridged):**
```json
{
  "isValid": true,
  "productSize": 106,
  "pairPenalty": 0.6983187292252637,
  "productTm": 74.01606258306629,
  "forward": { "sequence": "TCGTGAACTCGTGAACTCGT", "position": 76, "length": 20, "gcContent": 50.0, "meltingTemperature": 59.3 },
  "reverse": { "sequence": "CGGAACTTCGGAACTTCGGA", "position": 162, "length": 20, "gcContent": 55.0, "meltingTemperature": 60.0 }
}
```

With `"product_size_range": "150-250 100-149", "num_return": 3` the pairs end at 229, 230 and 237 (products 154/155/162); with `"pick_internal_oligo": true` primer3-py and this tool return PRIMER_RIGHT_0 = [197,20] and PRIMER_INTERNAL_0 = [154,24] `TCCGAAGTTCCGAAGTTCCGAAGT` (penalty 5.8754).

### Example 2: Invalid target region

`design_primers("ACGT…", 150, 100)` (start ≥ end) throws `ArgumentException`.

## See Also

- [generate_primer_candidates](generate_primer_candidates.md), [evaluate_primer](evaluate_primer.md), [primer_dimer](primer_dimer.md)
