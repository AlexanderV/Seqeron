---
type: source
title: "Evidence: EPIGEN-DMR-001 (Differentially methylated region detection)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-DMR-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-DMR-001-Evidence.md
source_commit: 547b4cdb6c5e5e4454d7cbdd9005a4741adc00e0
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-DMR-001

The validation-evidence artifact for test unit **EPIGEN-DMR-001** — **differentially methylated
region (DMR) detection** by the methylKit **tiling-window + Fisher's-exact-test** model. This is the
**fifth ingested unit of the Epigenetics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm is synthesized in its own
concept, [[differentially-methylated-regions]]; [[test-unit-registry]] tracks the unit. It compares
methylation between two samples, consuming the per-CpG β-values that
[[bisulfite-methylation-calling]] produces.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Akalin et al. (2012)** "methylKit", *Genome Biology* 13:R87 (PMC3491415, rank 1/3) — DMC/DMR
    definition (a base/region "classified as a differentially methylated cytosine (DMC) or region
    (DMR)" when the null is rejected); per-site methylation = **ratio of C/(C+T)** at each base;
    default cutoffs **q-value < 0.01 and %methylation difference > 25%**; hyper = higher methylation
    than control, hypo = lower.
  - **methylKit `tileMethylCounts`** man page (al2na/methylKit, rank 3) — defaults **`win.size=1000`,
    `step.size=1000`, `cov.bases=0`, `mc.cores=1`**; "summarizes methylated/unmethylated base counts
    over tiling windows across genome".
  - **methylKit `get.methylDiff`** source (diffMeth.R, rank 3) — signature
    `get.methylDiff(.Object, difference=25, qvalue=0.01, type="all")`; hyper
    `qvalue<qvalue & meth.diff > difference`, hypo `qvalue<qvalue & meth.diff < -1*difference`;
    **strict** greater-than, `difference` in percentage points.
  - **methylKit `calculateDiffMeth`** man page (rank 3) — one sample per group ⇒ **Fisher's exact
    test**; multiple samples per group ⇒ logistic regression (out of scope).
  - **Fisher's exact test** — Wikipedia citing Fisher 1922/1935 (rank 4, primary math) — single 2×2
    hypergeometric probability `p = (a+b)!(c+d)!(a+c)!(b+d)!/(a!b!c!d!n!)`; two-sided p = sum of
    probabilities of all same-margin tables with probability ≤ the observed table's; worked example
    a=1,b=9,c=11,d=3, n=24 → single-table **p ≈ 0.001346076**.

- **Documented corner cases / failure modes:** empty input → no tiles → **no DMRs**; sparse/low-
  coverage tiles filtered by `cov.bases` (default 0 keeps all; a tile needs ≥1 informative cytosine);
  one sample per group → Fisher (replicates → logistic regression, out of scope); **zero-coverage
  group** (row total 0) → degenerate 2×2, only the observed table feasible → two-sided **p = 1.0**;
  identical proportions / zero marginal → **p = 1.0**.

- **Datasets (documented oracles):**
  - *Fisher single-table probability* — a=1, b=9, c=11, d=3, n=24 → **≈0.001346076** (validates the
    hypergeometric term inside the p-value).
  - *Hyper-methylated window* — group1 (control) level 0.0 / coverage 20, group2 (treatment) level
    1.0 / coverage 20, 3 sites → pooled 2×2 methylated {g1=0, g2=60} / unmethylated {g1=60, g2=0};
    **meth.diff = +100**, two-sided Fisher **p ≈ 0** (complete separation); classified
    **Hypermethylated** (meth.diff > 25, q < 0.01).

- **Test-coverage recommendations:** MUST — hyper window (positive meth.diff + "Hypermethylated"),
  hypo window (negative meth.diff + "Hypomethylated"), |meth.diff| at/under cutoff NOT reported
  (strict `>`), tiling boundary (sites > win.size apart → separate windows), Fisher single-table
  probability matches ≈0.001346076, empty input → no DMRs. SHOULD — too-few-covered-sites not
  reported; degenerate group (zero coverage) → p 1.0, not reported. COULD — determinism (identical
  DMR list ordered by start).

## Deviations and assumptions

- **ASSUMPTION (per-window pooling):** methylKit's Fisher path operates on per-base counts; this unit
  pools the covered cytosines in a tile into one 2×2 table (Σ methylated vs Σ unmethylated per group)
  and applies Fisher's exact test to the window — mirrors `tileMethylCounts` (sum counts over the
  window) followed by Fisher's exact on the tiled counts, the documented methylKit tile→test pipeline
  (evidence-backed, not a free assumption).
- **ASSUMPTION (count reconstruction):** numC/numT are derived from `round(level × coverage)` and
  `round((1−level) × coverage)` because the repository `MethylationSite` stores a fractional level +
  integer coverage rather than raw C/T counts (methylKit stores numCs/numTs directly) — a
  representation detail, not a change to the C/(C+T) definition.

No source contradictions.
