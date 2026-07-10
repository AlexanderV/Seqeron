---
type: source
title: "Evidence: PRIMER-TM-001 (Primer3 weighted per-primer penalty objective function)"
tags: [validation, primer]
doc_path: docs/Evidence/PRIMER-TM-001-Evidence.md
sources:
  - docs/Evidence/PRIMER-TM-001-Evidence.md
source_commit: 92f89a5d97d45d2bac16b2d94467d1ee038879f7
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PRIMER-TM-001

The validation-evidence artifact for the **base** test unit **PRIMER-TM-001** — the
**Primer3 weighted per-primer penalty (objective function)** `p_obj_fn`, the score by which
Primer3 ranks candidate primers. Despite the unit ID, this file documents a **scoring /
selection** algorithm, not a melting-temperature calculation: Tm is only one input *term*.
The algorithm is synthesized in [[primer3-weighted-penalty-objective]]; the sibling
[[primer-tm-001-dimer-evidence]] (concept [[primer-dimer-thermodynamics-tm]]) records the
separate self-/hetero-dimer Tm extension carrying the same unit ID. This is one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; see
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Primer3 `libprimer3.cc` `p_obj_fn`** (rank 3, reference implementation) — the verbatim
    per-primer penalty: a **weighted sum of sign-gated one-sided deviations** accumulated term
    by term. Tm terms `temp_gt·(Tm−opt_tm)` / `temp_lt·(opt_tm−Tm)`; GC% terms
    `gc_content_gt·(GC−opt_gc)` / `gc_content_lt·(opt_gc−GC)`; size terms
    `length_gt·(len−opt_size)` / `length_lt·(opt_size−len)`; non-thermodynamic alignment
    `compl_any·self_any` + `compl_end·self_end`; `num_ns·<#N>`. **Lower is better** —
    `h.quality = p_obj_fn(...)`, primers sorted ascending. `gc_content` is a **percent**
    (`100·num_gc/num_gcat`). Defaults (`pr_set_default_global_args_2`): `opt_size=20`,
    `opt_tm=60.0`, `min_tm=57`, `max_tm=63`, `min_size=18`, `max_size=27`; default weights
    `temp_gt=temp_lt=length_gt=length_lt=1`, `gc_content_*=compl_*=num_ns=0`. So the
    **default** objective is `|Tm−60| + |len−20|`.
  - **Primer3 `libprimer3.h`** (rank 3) — `#define DEFAULT_OPT_GC_PERCENT PR_UNDEFINED_INT_OPT`;
    with GC weights 0 by default, GC% contributes nothing to the default objective (the manual's
    user-facing default `PRIMER_OPT_GC_PERCENT = 50.0` applies only when GC weights are set).
  - **Primer3 manual §19** "How Primer3 calculates the penalty value" (rank 2) — the
    user-facing formula (identical structure; right primers identical to left), "this penalty
    is the only score by which Primer3 evaluates the primers", published defaults
    (`PRIMER_OPT_GC_PERCENT=50.0`, `PRIMER_WT_SELF_ANY=0.0`, `PRIMER_MAX_SELF_ANY=8.00`), and
    that `SELF_ANY`/`SELF_END` are Primer3 `dpal` **local-alignment** scores (+1.00 complementary,
    −0.25 base·N, −1.00 mismatch, −2.00 gap; non-negative), **not** a closed form.
  - **Untergasser et al. 2012**, *Primer3 — new capabilities* (rank 1, Nucleic Acids Res
    40(15):e115) — "acceptable pairs [sorted] by a penalty function"; "'best' … operationally
    defined as minimizing a penalty function".
- **Datasets (hand-computed worked oracles):**
  - *Default weights* (`|Tm−60|+|len−20|`): (60,20,GC50,N0)→**0.0**; (63,20)→**3.0**;
    (57,18)→**5.0**; (62.5,22)→**4.5**.
  - *Non-default weights:* GC60 with `GC_GT=0.5`→**5.0**; GC40 with `GC_LT=0.5`→**5.0**;
    selfAny4 with `SELF_ANY=0.1`→**0.4**; selfEnd3 with `SELF_END=0.2`→**0.6**; N=2 with
    `NUM_NS=1`→**2.0**; combined (Tm62,len22,GC55,selfAny2,N1)→**8.0**.
- **Corner cases / failure modes:**
  - **Sign-gated terms:** at `Tm==opt` neither Tm term fires (0 contribution); same for size.
  - **Weight-zero short-circuit:** each `if(weights.X && …)` skips a 0-weight term (GC / self /
    num_ns default 0).
  - **Non-negative result:** every term is weight·(non-negative deviation) ⇒ penalty ≥ 0; 0 iff
    every parameter is at its optimum.
  - **Internal oligo** uses `o_args` weights (out of scope — PCR primers only, `p_args`).

## Deviations and assumptions

**ASSUMPTION:** `self_any` / `self_end` are **caller-supplied** alignment scores. The penalty
*arithmetic* (`WT_SELF_ANY·SELF_ANY`, etc.) is reproduced exactly, but the *value* of `SELF_ANY`
comes from Primer3's `dpal` local-alignment engine, which is not reproduced in this method. Under
**default** weights (`WT_SELF_ANY=WT_SELF_END=0`) these terms vanish, so the default objective is
fully reproduced; only the alignment-score *input* is external. No contradictions among sources.

Recommended coverage (MUST): default penalty `=|Tm−60|+|len−20|` (A–D); GC `_gt`/`_lt` fire only on
the correct side of `OPT_GC` (E,F); `SELF_ANY`/`SELF_END`/`NUM_NS` scale linearly with weight
(G,H,I); combined multi-term penalty (J); penalty at optimum exactly 0 and always ≥ 0. SHOULD:
`Tm==opt` and `len==opt` contribute 0 (sign gate); default-weights struct matches sourced constants.
COULD: lower penalty selects the better of two candidates (ordering semantics).
