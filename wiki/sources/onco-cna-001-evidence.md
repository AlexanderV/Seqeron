---
type: source
title: "Evidence: ONCO-CNA-001 (copy-number alteration classification — log2 ratio → absolute CN → CNA state)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CNA-001-Evidence.md
sources:
  - docs/Evidence/ONCO-CNA-001-Evidence.md
source_commit: f794f2454365bdf09e3440ac95c4f59db1ff0553
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CNA-001

The validation-evidence artifact for test unit **ONCO-CNA-001** — **copy-number alteration
classification**: a single **log2 copy ratio** → **absolute integer copy number** → discrete **CNA
state** (DeepDeletion / Loss / Neutral / Gain / Amplification). The **eighth ingested unit of the
Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[copy-number-alteration-classification]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (five, mutually consistent, no contradictions):**
  - **CNVkit `cnvlib/call.py`** (rank 3, reference implementation) — `_log2_ratio_to_absolute_pure`:
    `n = ref_copies · 2^log2`, diploid `ref_copies = 2` → **`n = 2·2^log2`**. `absolute_threshold`
    hard-threshold caller: CN = index of the **first** threshold with `log2 <= thresh`
    (**boundary-inclusive → lower bin**); above the last threshold `CN = ceil(2·2^log2)`. **Default
    thresholds `(−1.1, −0.25, 0.2, 0.7)` → states `[0,1,2,3,4+]`**. Docstring heuristic (50% clonality,
    ±0.1 noise): `DEL(0) < −1.1`, `LOSS(1) < −0.25`, `GAIN(3) ≥ +0.2`, `AMP(4) ≥ +0.7`. **NaN log2 →
    replaced with neutral CN** (no-call → diploid).
  - **CNVkit docs — `call` threshold method** (rank 3) — mapping table; cutoffs derive from
    `np.log2((copy_nums + .5) / 2)`; **germline-tuned `−0.4 / 0.3`** shown on the docs page vs the
    **tumor `−0.25 / 0.2`** source-code default the implementation follows; **purity ≥ 30% caveat**.
  - **GISTIC2.0 — Mermel et al. 2011, Genome Biology** (rank 1, peer-reviewed) — low-amplitude **noise
    threshold ±0.1 log2** (corroborates CNVkit's neutral bin) and high-amplitude cutoffs **+0.848 /
    −0.737** for amplification/deletion (consistent with AMP ≥ +0.7 / DEL ≤ −1.1).
  - **GISTIC2 docs — `-ta`/`-td`** (rank 2) — amplification/deletion thresholds default 0.1.
  - **Seqeron in-repo `StructuralVariantAnalyzer.DetectCNV`/`SegmentCopyNumber` (SV-CNV-001)** (rank 3,
    overlap check) — SV-CNV-001 does depth→log2→**integer** CN via `round(2·2^log2)` + segment merge but
    does **not** classify discrete CNA states; ONCO-CNA-001 is the complementary **oncology
    classification layer** (same conversion formula, CNVkit hard-threshold binning instead of `round`).

- **Documented corner cases / failure modes:** NaN log2 → neutral CN 2 (no-call); boundary inclusivity
  (`log2 <= thresh` → lower state, e.g. log2 = 0.7 → CN 3); Amplification is **open-ended**
  `ceil(2·2^log2)` (not fixed); GISTIC |log2| ≤ 0.1 noise band subsumed by the neutral bin.

- **Datasets (worked, deterministic — default thresholds −1.1, −0.25, 0.2, 0.7):**
  - **State/boundary table:** log2 −2.0 & −1.1 → CN 0 DeepDeletion; −1.0 & −0.25 → 1 Loss; 0.0 & 0.2 →
    2 Neutral; log2(3/2)=0.585 & 0.7 → 3 Gain; 1.0 → 4 & 2.0 → 8 Amplification.
  - **Absolute-CN formula `n = 2·2^log2`:** 0→2.0, 1→4.0, −1→1.0, log2(3/2)→3.0, 0.8→3.482…→ceil 4.

- **Coverage recommendations (7 items):** MUST-test each of the five CNA states for a representative
  in-bin log2; boundary inclusivity (−1.1→0, −0.25→1, 0.2→2, 0.7→3); the `n = 2·2^log2` formula; the
  above-last-threshold `ceil(2·2^log2)` amplification path (1.0→4, 2.0→8); NaN log2 → Neutral/CN 2.
  SHOULD-test batch order/length preservation. COULD-test custom-threshold override.

## Deviations and assumptions

- **ASSUMPTION — diploid (autosomal) reference ploidy = 2**, the CNVkit default (`ref_copies = ploidy =
  2`); sex-chromosome / haploid references out of scope. Source-backed default but not caller-overridable,
  so a non-diploid baseline would change every output.
- **Threshold-variant note (not a contradiction):** the CNVkit docs page's `−0.4 / 0.3` is the
  germline-precision alternative; the code default `−0.25 / 0.2` is the tumor-sample heuristic the
  implementation follows — the same docstring labels both.

No source contradictions — CNVkit (formula, binning, thresholds) and GISTIC2 (independent ±0.1 noise
band + high-amplitude cutoffs) corroborate one another across the whole classification.
