# Evidence Artifact: PRIMER-TM-001

**Test Unit ID:** PRIMER-TM-001
**Algorithm:** Primer3 weighted per-primer penalty (objective function)
**Date Collected:** 2026-06-23

---

## Online Sources

### Primer3 source — `libprimer3.cc`, `p_obj_fn` (the per-primer objective function)

**URL:** https://raw.githubusercontent.com/primer3-org/primer3/main/src/libprimer3.cc
**Accessed:** 2026-06-23
**Authority rank:** 3 (reference implementation — the canonical Primer3 source)

**Key Extracted Points:**

1. **Per-primer penalty is `p_obj_fn`:** For a left/right primer (`j == OT_LEFT || j == OT_RIGHT`) the penalty `sum` is accumulated term by term. Each term is added only when its weight is non-zero and the deviation has the right sign. Verbatim core terms:
   - `if (weights.temp_gt && h->temp > opt_tm) sum += weights.temp_gt * (h->temp - opt_tm);`
   - `if (weights.temp_lt && h->temp < opt_tm) sum += weights.temp_lt * (opt_tm - h->temp);`
   - `if (weights.gc_content_gt && h->gc_content > opt_gc_content) sum += weights.gc_content_gt * (h->gc_content - opt_gc_content);`
   - `if (weights.gc_content_lt && h->gc_content < opt_gc_content) sum += weights.gc_content_lt * (opt_gc_content - h->gc_content);`
   - `if (weights.length_lt && h->length < opt_size) sum += weights.length_lt * (opt_size - h->length);`
   - `if (weights.length_gt && h->length > opt_size) sum += weights.length_gt * (h->length - opt_size);`
   - (non-thermodynamic alignment) `if (weights.compl_any) sum += weights.compl_any * h->self_any;` and `if (weights.compl_end) sum += weights.compl_end * h->self_end;`
   - `if (weights.num_ns) sum += weights.num_ns * h->num_ns;`
2. **Lower is better:** `p_obj_fn` returns the non-negative `sum`; primers/pairs are sorted ascending by this value. Each acceptable primer's `h.quality = p_obj_fn(pa, &h, oligo->type)`.
3. **`gc_content` units are percent (0–100):** `h->gc_content = 100.0 * ((double)num_gc)/num_gcat;` (line 3856). So the GC term operates on percentages, like the Tm term operates on °C and the size term on bases.
4. **Default optima (`pr_set_default_global_args` / `_2`):** `opt_size = 20` (bases), `opt_tm = 60.0` (°C), `min_tm = 57.0`, `max_tm = 63.0`, `min_size = 18`, `max_size = 27`. `opt_gc_content = DEFAULT_OPT_GC_PERCENT`.
5. **Default weights (`p_args.weights`):** `temp_gt = 1`, `temp_lt = 1`, `length_gt = 1`, `length_lt = 1`, `gc_content_gt = 0`, `gc_content_lt = 0`, `compl_any = 0`, `compl_end = 0`, `num_ns = 0`, `pos_penalty = 1`, `seq_quality = 0`, `end_stability = 0`, `repeat_sim = 0`, `template_mispriming = 0`. Therefore the **default** objective for a primer with no N, no position penalty, is `temp_gt|temp_lt·|Tm−60| + length_gt|length_lt·|len−20|`.

### Primer3 source — `libprimer3.h` (`DEFAULT_OPT_GC_PERCENT`)

**URL:** https://raw.githubusercontent.com/primer3-org/primer3/main/src/libprimer3.h
**Accessed:** 2026-06-23
**Authority rank:** 3

**Key Extracted Points:**

1. **`#define DEFAULT_OPT_GC_PERCENT PR_UNDEFINED_INT_OPT`** — in the source, the default optimal GC% is *undefined*; combined with `gc_content_*` weights defaulting to 0, GC% contributes nothing to the default objective. The published manual documents the user-facing default of `PRIMER_OPT_GC_PERCENT = 50.0` (see below).

### Primer3 manual — §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE"

**URL:** https://raw.githubusercontent.com/primer3-org/primer3/main/src/primer3_manual.htm (rendered: https://primer3.org/manual.html)
**Accessed:** 2026-06-23
**Authority rank:** 2 (official specification/manual)

**Key Extracted Points:**

1. **Documented per-primer formula (verbatim structure):** `PRIMER_LEFT_4_PENALTY =` sum of, when `TM > OPT_TM`: `+ PRIMER_WT_TM_GT * (TM - OPT_TM)`; when `TM < OPT_TM`: `+ PRIMER_WT_TM_LT * (OPT_TM - TM)`; when `GC% > OPT_GC_PERCENT`: `+ PRIMER_WT_GC_PERCENT_GT * (GC% - OPT_GC_PERCENT)`; when `GC% < OPT_GC_PERCENT`: `+ PRIMER_WT_GC_PERCENT_LT * (OPT_GC_PERCENT - GC%)`; when `length > OPT_SIZE`: `+ PRIMER_WT_SIZE_GT * (length - OPT_SIZE)`; when `length < OPT_SIZE`: `+ PRIMER_WT_SIZE_LT * (OPT_SIZE - length)`; always `+ PRIMER_WT_SELF_ANY * SELF_ANY + PRIMER_WT_SELF_END * SELF_END + PRIMER_WT_NUM_NS * <number of N>`. "Right Primers (identical to Left Primers)."
2. **Selection meaning:** "In the second step, Primer3 calculates a penalty for each primer. This penalty is the only score by which Primer3 evaluates the primers." A lower number indicates a better primer/pair.
3. **Published default values:** `PRIMER_OPT_GC_PERCENT (float; default 50.0)`, `PRIMER_WT_SELF_ANY (float; default 0.0)`, `PRIMER_MAX_SELF_ANY (default 8.00)`.
4. **`SELF_ANY`/`SELF_END` are local-alignment scores:** "The scoring system gives 1.00 for complementary bases, -0.25 for a match of any base (or N) with an N, -1.00 for a mismatch, and -2.00 for a gap … Scores are non-negative." These are computed by Primer3's `dpal` local-alignment engine, not a closed-form formula.

### Untergasser et al. 2012, "Primer3—new capabilities and interfaces", Nucleic Acids Res 40(15):e115

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3424584/ (DOI https://doi.org/10.1093/nar/gks596)
**Accessed:** 2026-06-23
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. "Primer3 evaluates the primers and primer pairs according to various constraints and sorts acceptable pairs by a penalty function."
2. "the notion of 'best' primers or primer pairs is operationally defined as minimizing a penalty function." The paper defers the per-term details to the Primer3 manual.

---

## Documented Corner Cases and Failure Modes

### From `libprimer3.cc` `p_obj_fn`

1. **Sign-gated terms:** A term is only added when the deviation has the matching sign; e.g. when `Tm == opt_tm` neither the `temp_gt` nor `temp_lt` term fires, so a primer exactly at optimum contributes 0 from that parameter.
2. **Weight-zero short-circuit:** Each `if (weights.X && ...)` short-circuits when the weight is 0, so a 0-weight term never contributes (matters because GC/self/num_ns default to 0).
3. **Non-negative result:** Every term is weight·(non-negative deviation), so the total penalty is ≥ 0; 0 means every parameter is at its optimum.

### From the manual §19

1. **Internal vs primer:** Left and right primers use identical formulas (`p_args`); the internal oligo uses `o_args` weights — out of scope here (PCR primers only).

---

## Test Datasets

### Dataset: Hand-computed worked examples from the published formula (default weights/optima)

**Source:** Primer3 `libprimer3.cc` `p_obj_fn` + default weights/optima (Source above). Default primer weights: `WT_TM_GT=WT_TM_LT=WT_SIZE_GT=WT_SIZE_LT=1`, `WT_GC=WT_SELF=WT_NUM_NS=0`; optima `OPT_TM=60`, `OPT_SIZE=20`, `OPT_GC=50`.

| # | Inputs (Tm°C, len, GC%, N) | Default penalty = |Tm−60|·1 + |len−20|·1 | Value |
|---|----------------------------|------------------------------------------|-------|
| A | Tm=60, len=20, GC=50, N=0 | 0 + 0 | 0.0 |
| B | Tm=63, len=20 | 3·1 + 0 | 3.0 |
| C | Tm=57, len=18 | 3·1 + 2·1 | 5.0 |
| D | Tm=62.5, len=22 | 2.5·1 + 2·1 | 4.5 |

### Dataset: Non-default weights (exercise GC, self, num_ns terms)

**Source:** same formula, weights set explicitly per the manual's parameters.

| # | Inputs | Weights | Penalty | Value |
|---|--------|---------|---------|-------|
| E | Tm=60, len=20, GC=60 | GC_GT=0.5 | 0.5·(60−50) | 5.0 |
| F | Tm=60, len=20, GC=40 | GC_LT=0.5 | 0.5·(50−40) | 5.0 |
| G | Tm=60, len=20, selfAny=4 | SELF_ANY=0.1 | 0.1·4 | 0.4 |
| H | Tm=60, len=20, selfEnd=3 | SELF_END=0.2 | 0.2·3 | 0.6 |
| I | Tm=60, len=20, N=2 | NUM_NS=1 | 1·2 | 2.0 |
| J | Tm=62, len=22, GC=55, selfAny=2, N=1 | TM_GT=1, SIZE_GT=1, GC_GT=0.5, SELF_ANY=0.25, NUM_NS=1 | 1·2 + 1·2 + 0.5·5 + 0.25·2 + 1·1 | 8.0 |

---

## Assumptions

1. **ASSUMPTION: `self_any` / `self_end` alignment scores are caller-supplied to the penalty method.** The penalty term `WT_SELF_ANY * SELF_ANY` is reproduced exactly, but the *value* of `SELF_ANY` comes from Primer3's `dpal` local-alignment engine (scoring +1.00 complementary / −1.00 mismatch / −2.00 gap), which is not reproduced here. Under the **default** weights (`WT_SELF_ANY = WT_SELF_END = 0`) these terms contribute nothing, so the default objective is fully reproduced; for non-default weights the caller supplies the alignment score. This is the only term whose *input* is not computed in-method; the *penalty arithmetic* is exact.

---

## Recommendations for Test Coverage

1. **MUST Test:** Default-weight penalty equals `|Tm−60| + |len−20|` for datasets A–D — Evidence: `libprimer3.cc` `p_obj_fn` + default weights.
2. **MUST Test:** GC `_gt`/`_lt` terms fire only on the correct side of `OPT_GC` (datasets E, F) — Evidence: `p_obj_fn` gc_content branches.
3. **MUST Test:** `SELF_ANY`, `SELF_END`, `NUM_NS` terms scale linearly with their weight (datasets G, H, I) — Evidence: `p_obj_fn` compl_any/compl_end/num_ns lines.
4. **MUST Test:** Combined multi-term penalty (dataset J) — Evidence: `p_obj_fn` (sum of terms).
5. **MUST Test:** Penalty at optimum is exactly 0; penalty is always ≥ 0 — Evidence: sign-gated non-negative terms.
6. **SHOULD Test:** `Tm == opt` and `len == opt` contribute 0 (sign gate) — Rationale: boundary documented in `p_obj_fn`.
7. **SHOULD Test:** Default-weights struct matches the source defaults (TM/SIZE=1, GC/SELF/NUM_NS=0; OPT_TM=60, OPT_SIZE=20, OPT_GC=50) — Rationale: guards against drift of sourced constants.
8. **COULD Test:** Lower penalty selects the better of two candidates (ordering) — Rationale: "lower is better" selection semantics.

---

## References

1. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research 40(15):e115. https://doi.org/10.1093/nar/gks596
2. Koressaar T, Remm M (2007). Enhancements and modifications of primer design program Primer3. Bioinformatics 23(10):1289–1291. https://doi.org/10.1093/bioinformatics/btm091
3. Primer3 source code, `src/libprimer3.cc` (`p_obj_fn`, `pr_set_default_global_args_2`) and `src/libprimer3.h` (`DEFAULT_OPT_GC_PERCENT`), branch `main`. https://github.com/primer3-org/primer3
4. Primer3 manual, §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE" and global input tags (default weights/optima). https://primer3.org/manual.html

---

## Change History

- **2026-06-23**: Initial documentation (PRIMER-TM-001 — Primer3 weighted penalty objective).
