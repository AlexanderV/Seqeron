---
type: source
title: "Evidence: SV-CNV-001 (read-depth copy-number variation detection)"
tags: [validation, structural-variant]
doc_path: docs/Evidence/SV-CNV-001-Evidence.md
sources:
  - docs/Evidence/SV-CNV-001-Evidence.md
source_commit: 59811dacff3428aa9f6ae78b68795bba34ce864d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SV-CNV-001

The validation-evidence artifact for test unit **SV-CNV-001** — **Read-Depth Copy Number Variation
Detection** (windowed read depth → log2 ratio → integer copy number → deletion/duplication call). The
**read-depth-segmentation member of the germline structural-variant (SV) family** and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct
method is synthesized in [[read-depth-cnv-segmentation]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the read-depth CNV paradigm):**
  - **Yoon, Xuan, Makarov, Ye & Sebat 2009** (Genome Research 19(9):1586–1592, rank 1 primary paper;
    open-access PMC2752127) — the **RD ∝ CN** core hypothesis ("a linear relationship between coverage
    and copy number" across CN 1/2/3); **windowed counting** ("counting the number of mapped reads in
    100-bp windows, assigning each read only once by its start position"); and the **GC-correction
    formula** verbatim `r_i' = r_i · m / m_GC` (per-window counts normalised to the overall median `m`).
  - **CNVkit `call.py`** (etal/cnvkit master, rank 3; primary Talevich et al. 2016 PLoS Comput Biol) —
    `_log2_ratio_to_absolute_pure` returns `n = ref_copies · 2^log2` then `round().astype(int)`, so for
    a diploid reference **CN = round(2·2^log2)**; `log2_ratio = log2(ncopies / ploidy)`; the anchors
    **log2(1/2) = −1.0** (single-copy loss), **log2(3/2) = 0.585** (single-copy gain), **log2(2/2) = 0**;
    the "threshold" method (`absolute_threshold`, default cutoffs `(−1.1, −0.25, 0.2, 0.7)`); and the
    non-negativity clamp `ncopies = max(0.0, ncopies)`.
  - **CNVkit "Calling copy number gains and losses" docs** (rank 3) — the diploid conversion "absolute
    copy number = 2 · 2^(log2 value)" (e.g. log2 0.38 → 2·2^0.38 = 2.6).

- **Documented corner cases / failure modes:**
  - **Zero / very-low-coverage window** — RD = 0 gives `log2(0/ref)` undefined (−∞); treated as an
    **unobserved / no-call** window (homozygous-deletion candidate), not a −∞ ratio.
  - **NaN log2** — CNVkit replaces a `nan` log2 with the **neutral** copy number (CN 2), a no-call rather
    than a spurious call.
  - **Negative extrapolation clamped** — CN is physically **≥ 0** (`max(0.0, ncopies)`).

- **Datasets (deterministic, derived from the cited formulas — reference/diploid baseline RD = 100,
  window size 4 positions):**
  - `100,100,100,100` → mean 100 → log2(1.0) = 0 → CN 2 → **Neutral**.
  - `50,50,50,50` → mean 50 → log2(0.5) = −1.0 → CN 1 → **Deletion (loss)**.
  - `0,0,0,0` → mean 0 → log2(0) undefined → **no-call** (homozygous-deletion candidate).
  - `150,150,150,150` → mean 150 → log2(1.5) = 0.585 → CN 3 → **Duplication (gain)**.
  - `200,200,200,200` → mean 200 → log2(2.0) = 1.0 → CN 4 → **Duplication (amplification)**.
  - CN check: round(2·2^−1) = 1, round(2·2^0.585) = 3, round(2·2^1) = 4, round(2·2^0) = 2.

- **Coverage recommendations (9 items):** MUST — RD = reference → log2 0 → CN 2 (Neutral); half
  reference → log2 −1 → CN 1 (loss); 1.5× → log2 0.585 → CN 3 (gain); 2× → log2 1.0 → CN 4; windowing
  summarises depth into non-overlapping windows with a per-window mean; a zero-depth window → no-call
  (not −∞). SHOULD — reference RD defaults to the overall median of window means when no baseline is
  given; CN is non-negative and the log2/CN relationship is monotonic non-decreasing in RD. COULD —
  empty input → empty output, null input → throws `ArgumentNullException`.

## Deviations and assumptions

- **ASSUMPTION — reference (diploid baseline) depth = the overall median of the windowed read depths.**
  Yoon normalises to the overall median `m`; CNVkit's log2 is taken against a reference profile. With no
  external reference supplied to `DetectCNV`, the per-sample reference RD is the overall median of the
  non-zero window means. Correctness-affecting (it sets the log2 = 0 anchor) but the source-supported
  self-reference choice; an explicit baseline can override it.
- **ASSUMPTION — diploid ploidy (2) is the copy-number baseline.** `CN = ploidy · 2^log2` and the anchors
  log2(1/2) = −1 / log2(3/2) = 0.585 are stated for a diploid genome; the standard human autosomal
  baseline used by both sources.

No source contradictions — Yoon 2009 (RD∝CN, windowed counting, GC correction) and CNVkit (log2 → CN
arithmetic, anchors, clamping) cover disjoint stages of one pipeline and are mutually consistent.
