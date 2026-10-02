# design_primers

Design a forward/reverse PCR primer pair flanking a target region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_primers` |
| **Method ID** | `PrimerDesigner.DesignPrimers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Implements Primer3's pair search (`libprimer3.cc` `choose_pair_or_triple` / `characterize_pair` / `obj_fn`; verified against primer3-py 2.3.1 `design_primers`):

- **Candidates.** Forward candidates end at or before `target_start`, reverse candidates (evaluated on the reverse complement) start at or after `target_end`; only positions that can give a product in `product_size_range` are enumerated. A candidate is kept when it passes every per-primer limit of `PrimerParameters` (length, GC %, Primer3 SantaLucia 1998 NN Tm at `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc` — Primer3 defaults 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM; the same conditions drive the ntthal structure / pair complementarity values and the product Tm, as in Primer3 — poly-X, dinucleotide repeat, no non-ACGT base).
- **Real filters of a pair** (in Primer3's order): product size in the current `product_size_range` range (default **100–300 bp**, Primer3 PRIMER_PRODUCT_SIZE_RANGE; ranges are tried in order), `|Tm_f − Tm_r| ≤ max_tm_difference` (default **5 °C**; Primer3's own default is 100), and — with the default `PrimerStructureScreen.Primer3Thermodynamic` — Primer3's thermodynamic structure screen: each primer's ntthal self-dimer (SELF_ANY_TH), 3′ self-dimer (SELF_END_TH) and hairpin (HAIRPIN_TH) Tm ≤ **47 °C**, and the pair's ntthal hetero-dimer (COMPL_ANY_TH) and 3′ hetero-dimer (COMPL_END_TH) Tm ≤ **47 °C** (`MaxStructureTm`). With `StructureScreen = Heuristic` the screen is the sequence-only stem-loop test plus the alignment 3′ complementarity test. With `pick_internal_oligo` the pair also needs an internal hybridization oligo strictly between the primers (Primer3 PRIMER_INTERNAL_* defaults).
- **Selection.** The pair with the lowest Primer3 pair penalty (PRIMER_PAIR_WT_PR_PENALTY = 1: sum of the per-primer penalties |Tm − OptimalTm| + |length − OptimalLength|); ties as Primer3's `compare_primer_pair`. `pairs` lists up to `num_return` ranked pairs (PRIMER_NUM_RETURN; each later pair excludes the earlier ones). When no pair qualifies, the individually best primers are returned with `isValid = false` and a message naming the violated constraint. `product_size = reverse.Position + reverse.Length − forward.Position`.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L104)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | DNA template (A/C/G/T), non-empty. |
| `target_start` | integer | Yes | 0-based inclusive start of the target region (≥ 0). |
| `target_end` | integer | Yes | 0-based **exclusive** end of the target region (`target_start < target_end < template.Length`); primers never overlap `[target_start, target_end)`. |
| `parameters` | object | No | Optional `PrimerParameters`; defaults (18–25 nt, 40–60% GC, 57–63 °C Tm, poly-X 4, dinucleotide repeat 4, Primer3 thermodynamic structure screen at 47 °C) when null. |
| `product_size_range` | string | No | PRIMER_PRODUCT_SIZE_RANGE in Primer3 syntax, e.g. `"100-300"` or `"150-250 100-400"` (default `100-300`). |
| `max_tm_difference` | number | No | PRIMER_PAIR_MAX_DIFF_TM in °C (default 5). |
| `num_return` | integer | No | PRIMER_NUM_RETURN: number of ranked pairs in `pairs` (default 1, ≥ 1). |
| `pick_internal_oligo` | boolean | No | PRIMER_PICK_INTERNAL_OLIGO (default false). |
| `pair_max_compl_any` / `pair_max_compl_end` | number | No | PRIMER_PAIR_MAX_COMPL_ANY / _COMPL_END (default 8 / 3); used only with `parameters.StructureScreen = Primer3Alignment` (Primer3 PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0: per-primer dpal self_any ≤ `parameters.MaxSelfAny` 8, self_end ≤ `MaxSelfEnd` 3; internal oligo ≤ 12 / 12). |
| `salt_monovalent` / `salt_divalent` / `dntp_conc` / `dna_conc` | number | No | PRIMER_SALT_MONOVALENT (mM, > 0; default 50) / PRIMER_SALT_DIVALENT (mM, ≥ 0; 1.5) / PRIMER_DNTP_CONC (mM, ≥ 0; 0.6) / PRIMER_DNA_CONC (nM, > 0; 50); override `parameters`. Illegal values → `ArgumentOutOfRangeException` (Primer3 `_pr_data_control`). |
| `opt_gc_percent` / `wt_gc_percent_gt` / `wt_gc_percent_lt` | number | No | PRIMER_OPT_GC_PERCENT (default 50) and PRIMER_WT_GC_PERCENT_GT / _LT (default 0) of the per-primer penalty. |
| `internal_salt_monovalent` / `internal_salt_divalent` / `internal_dntp_conc` / `internal_dna_conc` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT / _SALT_DIVALENT / _DNTP_CONC / _DNA_CONC of the internal oligo (defaults 50 mM / 0 / 0 / 50 nM). |
| `internal_opt_gc_percent` / `internal_wt_gc_percent_gt` / `internal_wt_gc_percent_lt` | number | No | PRIMER_INTERNAL_OPT_GC_PERCENT (default 50) and PRIMER_INTERNAL_WT_GC_PERCENT_GT / _LT (default 0). Verified (all condition arguments) against primer3-py 2.3.1 `design_primers`: 2500/2500 random templates. |
| `gc_clamp` / `max_end_gc` / `max_end_stability` | integer / integer / number | No | Primer3 3′-end checks of each primer: PRIMER_GC_CLAMP (consecutive 3′ G/C required; default 0; ≤ minimum primer length), PRIMER_MAX_END_GC (max G/C in the last 5 bases, 0–5; default 5), PRIMER_MAX_END_STABILITY (max `end_stability` = −ΔG of the 3′ pentamer, kcal/mol, ≥ 0; default 100); override `parameters`. Illegal values → `ArgumentOutOfRangeException`. |
| `min_left_three_prime_distance` / `min_right_three_prime_distance` / `min_three_prime_distance` | integer | No | PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE (≥ −1, default −1 = primers may be reused across the `pairs`): after a pair is selected, later pairs may not use a left/right primer whose 3′ end is fewer bases away (0 = not the identical primer); `min_three_prime_distance` (PRIMER_MIN_THREE_PRIME_DISTANCE) sets both and cannot be combined with the specific ones. Verified (with the 3′-end checks) against primer3-py 2.3.1 `design_primers`: 2200/2200 random templates (9590 pairs). |
| `mispriming_library` | object | No | PRIMER_MISPRIMING_LIBRARY as a name → sequence object (primer3-py `misprime_lib`), e.g. `{"Alu*2": "GGCCGGGCGCGG…"}`. A `*weight` (0–100, C `strtod`) after the name scales that entry; sequences are upper-cased, whitespace removed, IUPAC codes kept, other characters → N. Each primer is aligned (Primer3 dpal: +1/−1/−0.25 N/−2 gap, max gap 1, anchored at the primer's 3′ end) with every entry and its reverse complement ("reverse <name>"); the weighted maximum is reported per primer and the pair score is the maximum over entries of the integer part of left + right. Empty or illegal-weight entries → `ArgumentException`. |
| `max_library_mispriming` / `pair_max_library_mispriming` | number | No | PRIMER_MAX_LIBRARY_MISPRIMING (default 12; compared as a C `short`, so 12.9 acts as 12) / PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (default 24): a primer / pair whose library score is strictly greater is rejected. |
| `wt_library_mispriming` / `pair_wt_library_mispriming` | number | No | PRIMER_WT_LIBRARY_MISPRIMING / PRIMER_PAIR_WT_LIBRARY_MISPRIMING (default 0): weight of the library score in the per-primer / pair penalty; non-zero without `mispriming_library` → `ArgumentException` (Primer3 `_pr_data_control`). |
| `lib_ambiguity_codes_consensus` | boolean | No | PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS (default false = Primer3 0): false — IUPAC codes in the library never match, N scores −0.25, and right primers are aligned without the 3′ anchor (Primer3 quirk); true — an IUPAC code matches every base it represents (runs of N match any primer). Verified against primer3-py 2.3.1 `design_primers(misprime_lib=…)`: see `PRIMER-DESIGN-001` F41 in `docs/Validation/review-2026-09/B07.md`. Also applies to `internal_mishyb_library` (Primer3 has one global setting). |
| `internal_mishyb_library` | object | No | PRIMER_INTERNAL_MISHYB_LIBRARY (with `pick_internal_oligo`): name → sequence object (primer3-py `mishyb_lib`, same format as `mispriming_library`); the internal oligo is aligned with every entry and its reverse complement by Primer3's dpal unanchored local alignment; reported as `internalOligo.libraryMishyb` / `libraryMishybName`. |
| `internal_max_library_mishyb` / `internal_wt_library_mishyb` | number | No | PRIMER_INTERNAL_MAX_LIBRARY_MISHYB (default 12; an internal oligo with a strictly greater score is rejected) / PRIMER_INTERNAL_WT_LIBRARY_MISHYB (default 0: weight of the score in the internal-oligo penalty, which enters the pair penalty through PRIMER_PAIR_WT_IO_PENALTY; non-zero without `internal_mishyb_library` → `ArgumentException`). Verified against primer3-py 2.3.1 `design_primers(mishyb_lib=…)`: see F42 in `docs/Validation/review-2026-09/B07.md`. |
| `thermodynamic_template_alignment` | boolean | No | PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT (default false = Primer3 0): false — template mispriming is Primer3's dpal score (3′-anchored local alignment of the primer with its own template strand 5′ and 3′ of its site, and with the opposite strand); true — the ntthal 3′-end (END1) Tm in °C at the primer conditions (template ≤ 10000 nt, else `ArgumentException`). Only the active mode's limits / weights apply. |
| `max_template_mispriming` / `max_template_mispriming_th` | number | No | PRIMER_MAX_TEMPLATE_MISPRIMING (alignment mode) / _TH (thermodynamic mode, °C); default −100 = not checked; a primer whose template mispriming value is strictly greater is rejected. > 32767 in alignment mode → `ArgumentException`. |
| `pair_max_template_mispriming` / `pair_max_template_mispriming_th` | number | No | PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING / _TH (default −100): pair value = max(left same-strand + right other-strand, left other-strand + right same-strand). Alignment mode checks a limit ≥ 0; thermodynamic mode, as Primer3, any non-zero limit — so `pair_wt_template_mispriming_th` > 0 with the default −100 rejects every pair; use 0 for "no limit". |
| `wt_template_mispriming` / `wt_template_mispriming_th` / `pair_wt_template_mispriming` / `pair_wt_template_mispriming_th` | number | No | PRIMER_[PAIR_]WT_TEMPLATE_MISPRIMING[_TH] (default 0, must be ≥ 0): weight in the per-primer / pair penalty (alignment: linear; thermodynamic: Primer3's 5 °C `temp_cutoff` rule). Verified against primer3-py 2.3.1 `design_primers`: see F43 in `docs/Validation/review-2026-09/B07.md`. |

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
| `internalOligo` | object \| null | PRIMER_INTERNAL_0_* (`sequence`, `start`, `length`, `tm`, `gcPercent`, `selfAnyTh`, `selfEndTh`, `hairpinTh`, `penalty`). |
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
