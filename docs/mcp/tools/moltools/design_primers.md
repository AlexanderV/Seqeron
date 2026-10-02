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

- **Candidates.** Forward candidates end at or before `target_start`, reverse candidates (evaluated on the reverse complement) start at or after `target_end`; only positions that can give a product in `product_size_range` are enumerated. A candidate is kept when it passes every per-primer limit of `PrimerParameters` (length, GC %, Primer3-default SantaLucia 1998 NN Tm at 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM, poly-X, dinucleotide repeat, no non-ACGT base).
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
