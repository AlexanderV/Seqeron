---
type: concept
title: "Primer3 weighted per-primer penalty (objective function)"
tags: [primer, algorithm, validation]
sources:
  - docs/algorithms/MolTools/Primer3_Penalty_Objective.md
  - docs/Evidence/PRIMER-TM-001-Evidence.md
source_commit: 540ba0d2f76748172d23905afb56af039fd7c75f
created: 2026-07-10
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: primer-tm-001-evidence
      evidence: "Test Unit ID: PRIMER-TM-001 ... Algorithm: Primer3 weighted per-primer penalty (objective function)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: primer-tm-001-evidence
      evidence: "Both carry unit ID PRIMER-TM-001; the penalty's temp_gt/temp_lt terms consume a per-primer Tm, and self_any/self_end consume dimer/self-alignment scores"
      confidence: high
      status: current
---

# Primer3 weighted per-primer penalty (objective function)

The **score by which Primer3 ranks candidate PCR primers** — its `p_obj_fn` **objective /
penalty function** (test unit **PRIMER-TM-001**, base unit). Despite the "TM" in the unit ID
this is a **selection / scoring** algorithm, not a melting-temperature calculation: Tm enters
only as one input *term*. It is the sibling of [[primer-dimer-thermodynamics-tm]] (the
self-/hetero-dimer Tm extension sharing the same unit ID) in the PCR primer-design
**PRIMER-\* / MolTools** family. The literature-traced record is [[primer-tm-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
artifact pattern.

## The penalty as a weighted sum of one-sided deviations

Primer3's two-step selection first filters by hard min/max constraints, then scores each
surviving primer by a **penalty** — the only value it sorts on (**lower is better**;
`h.quality = p_obj_fn(...)`, primers sorted ascending). The penalty is a **weighted sum of
sign-gated, one-sided deviations** from each parameter's optimum:

```
penalty =  temp_gt·(Tm − opt_tm)         when Tm > opt_tm
         + temp_lt·(opt_tm − Tm)         when Tm < opt_tm
         + gc_content_gt·(GC% − opt_gc)  when GC% > opt_gc
         + gc_content_lt·(opt_gc − GC%)  when GC% < opt_gc
         + length_gt·(len − opt_size)    when len > opt_size
         + length_lt·(opt_size − len)    when len < opt_size
         + compl_any·self_any + compl_end·self_end
         + num_ns·(#N bases)
```

Left and right primers use the **identical** formula and the same `p_args` weights (the
internal oligo uses a separate `o_args` set — out of scope here). `gc_content` is a **percent**
(`100·num_gc/num_gcat`), so the GC term operates on 0–100 like Tm operates on °C and size on
bases.

## Defaults reduce it to `|Tm−60| + |len−20|`

From `pr_set_default_global_args_2`: `opt_size=20`, `opt_tm=60.0` (with `min_tm=57`,
`max_tm=63`, `min_size=18`, `max_size=27`). Default weights are `temp_gt = temp_lt =
length_gt = length_lt = 1` and `gc_content_* = compl_any = compl_end = num_ns = 0`. So the
**default** objective collapses to

```
penalty = |Tm − 60| + |len − 20|
```

The optimal GC% is itself undefined in the header
(`#define DEFAULT_OPT_GC_PERCENT PR_UNDEFINED_INT_OPT`); with GC weights 0, GC% contributes
nothing by default. The published user-facing default `PRIMER_OPT_GC_PERCENT = 50.0` only takes
effect once a GC weight is set.

## Contract and structural properties

- **Sign-gated:** a term fires only when the deviation has the matching sign. At `Tm == opt_tm`
  (or `len == opt_size`) **neither** side's term fires — a primer exactly at optimum contributes
  0 from that parameter.
- **Weight-zero short-circuit:** each `if (weights.X && …)` skips a 0-weight term entirely, which
  is why GC / self-complementarity / num_ns cost nothing under defaults.
- **Non-negative:** every term is weight·(non-negative deviation) ⇒ penalty ≥ 0, and **penalty
  = 0 iff every parameter sits exactly at its optimum**.
- **`self_any` / `self_end` are inputs, not computed here:** these are Primer3 `dpal`
  local-alignment scores (+1.00 complementary, −0.25 base·N, −1.00 mismatch, −2.00 gap;
  non-negative). The penalty reproduces the *arithmetic* `WT_SELF_ANY·SELF_ANY` exactly but the
  score *value* is caller-supplied. Under default weights (both 0) this term vanishes, so the
  default objective is fully self-contained.

## Worked oracles

Default weights (`|Tm−60| + |len−20|`): (Tm 60, len 20) → **0.0**; (63, 20) → **3.0**;
(57, 18) → **5.0**; (62.5, 22) → **4.5**. Non-default weights: GC 60 with `GC_GT=0.5` → **5.0**;
selfAny 4 with `SELF_ANY=0.1` → **0.4**; N=2 with `NUM_NS=1` → **2.0**; combined
(Tm 62, len 22, GC 55, selfAny 2, N 1) with `TM_GT=SIZE_GT=NUM_NS=1, GC_GT=0.5, SELF_ANY=0.25`
→ `2 + 2 + 2.5 + 0.5 + 1` = **8.0**.

## Scope: which Primer3 terms this reproduces

The algorithm spec (`docs/algorithms/MolTools/Primer3_Penalty_Objective.md`,
`PrimerDesigner.CalculatePrimer3Penalty`) implements the **left/right-primer branch** of
`p_obj_fn` for the six core terms only (Tm, size, GC%, self_any, self_end, num_ns), as a
**deterministic O(1) weighted sum** — exact, not heuristic or probabilistic. Everything else
in Primer3's full objective is deliberately **out of scope**, and every omitted term carries a
default weight of 0, so the default penalty is unaffected:

- **Not implemented (per-primer):** the thermodynamic-alignment branch (`*_TH` terms with
  `temp_cutoff`), `pos_penalty`, `end_stability`, `seq_quality`, `repeat_sim`,
  `template_mispriming`. Rely on Primer3 itself if any of these weights are set non-zero.
- **Not implemented (pair-level):** the whole `PRIMER_PAIR_*` objective — Tm difference,
  product-size deviation, pair complementarity. This unit scores a single oligo; pair ranking
  is a separate (future) concern — currently handled by the greedy
  [[primer-design]] pair-selection pipeline (PRIMER-DESIGN-001), which still ranks candidates
  with the legacy additive `CalculatePrimerScore`, not this objective.

Because these all default to weight 0, the implementation still matches Primer3's *default*
selection exactly; the gap only appears when a caller enables an unsupported weight.

The unit also leaves the legacy convenience `Score` (`CalculatePrimerScore`, used by
`EvaluatePrimer` / `DesignPrimers`) **unchanged** alongside the new method — the Primer3-anchored
`CalculatePrimer3Penalty` is the validated objective, while `Score` remains for backward
compatibility. The Tm term input comes from the SEQ-THERMO-001 routine described in
[[melting-temperature]].
