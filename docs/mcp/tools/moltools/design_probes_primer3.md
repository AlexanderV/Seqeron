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

- Source: [ProbeDesigner.cs#L1121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1121)
- Algorithm doc: [Hybridization_Probe_Design.md](../../../algorithms/MolTools/Hybridization_Probe_Design.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | Template DNA (non-empty; case-insensitive; probes are picked on this strand). |
| `num_return` | integer | No | PRIMER_NUM_RETURN (≥ 0, default 5). |
| `min_size` / `opt_size` / `max_size` | integer | No | PRIMER_INTERNAL_MIN/OPT/MAX_SIZE (18 / 20 / 27; 1 ≤ min ≤ max ≤ 36). |
| `min_tm` / `opt_tm` / `max_tm` | number | No | PRIMER_INTERNAL_MIN/OPT/MAX_TM in °C (57 / 60 / 63). |
| `min_gc_percent` / `max_gc_percent` | number | No | PRIMER_INTERNAL_MIN/MAX_GC (20 / 80 %). |
| `max_poly_x` | integer | No | PRIMER_INTERNAL_MAX_POLY_X (5). |
| `max_self_any_th` / `max_self_end_th` / `max_hairpin_th` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH (47 °C). |
| `monovalent_mm` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT (50 mM). |
| `divalent_mm` | number | No | PRIMER_INTERNAL_SALT_DIVALENT (0 mM). |
| `dntp_mm` | number | No | PRIMER_INTERNAL_DNTP_CONC (0 mM). |
| `dna_conc_nm` | number | No | PRIMER_INTERNAL_DNA_CONC (50 nM). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `probes` | object[] | Best first. Each: `sequence`, `start` (0-based), `length`, `tm`, `gcPercent`, `selfAnyTh`, `selfEndTh`, `hairpinTh`, `penalty` (the Primer3 `PRIMER_INTERNAL_n_*` values). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Template sequence cannot be null or empty |
| 1002 | num_return cannot be negative |
| 1003 | Sizes must satisfy 1 ≤ min_size ≤ max_size ≤ 36 |

## Examples

### Example 1: Primer3 defaults

120-nt template `CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC`.

**Response (abridged):** first probe `start = 27`, `length = 22`, `tm = 57.5737`, `penalty = 4.4263`, `selfAnyTh = 15.4317`, `selfEndTh = 1.5232`, `hairpinTh = 37.4725` (= primer3-py `PRIMER_INTERNAL_0_*`); positions of the 5 probes (27,22), (27,23), (67,24), (68,23), (67,23).

## See Also

- [design_probes](design_probes.md), [validate_probe](validate_probe.md), [design_primers](design_primers.md)
