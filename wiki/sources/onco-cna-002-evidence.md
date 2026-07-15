---
type: source
title: "Evidence: ONCO-CNA-002 (focal amplification detection — GISTIC2 length-based focal/broad split + oncogene mapping)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CNA-002-Evidence.md
sources:
  - docs/Evidence/ONCO-CNA-002-Evidence.md
source_commit: c4b1520f4036a9d9cc96f01c35e49380b4dfa873
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CNA-002

The validation-evidence artifact for test unit **ONCO-CNA-002** — **focal amplification detection**:
keep segments that are both **amplified** (log2 gain > `t_amp` 0.1) and **focal** (segment length
< 98% of its chromosome arm), then **map** each focal amplification's arm to a panel of known
**oncogenes**. The **ninth ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in its own concept, [[focal-amplification-detection]]; [[test-unit-registry]] tracks the
unit.

## What this file records

- **Online sources (four, mutually consistent, no contradictions):**
  - **GISTIC2.0 — Mermel et al. 2011, Genome Biology 12:R41** (rank 1, peer-reviewed) — the
    **length-based focal/arm-level split**: focal SCNAs are those with "length < 98% of a chromosome
    arm"; events "occupying more than 98% of a chromosome arm" are arm-level; the procedure removes
    all such arm-level events "leaving only the focal events." Length "provides a natural basis for
    classifying events as 'arm-level' and 'focal' based purely on length" — GISTIC2.0's shift **away
    from amplitude-based** toward **length-based** filtering. (Historical amplitude context: ±0.1 log2
    eliminates only low-level artifacts; older 0.848/−0.737 high-amplitude cutoffs.)
  - **Broad Institute — GISTIC2 docs** (rank 3, parameter reference) — verbatim `broad_len_cutoff`
    default **0.98** ("threshold used to distinguish broad from focal events, given in units of
    fraction of chromosome arm"); `t_amp` default **0.1** ("regions with a copy number gain above this
    positive value are considered amplified"); `t_del` 0.1 recorded for completeness (deletions →
    ONCO-CNA-003, out of scope).
  - **CNVkit — calling docs** (rank 3) — a single-copy gain in a pure sample has copy ratio 3/2 →
    **log2(3/2) = 0.585**; any gain has log2 > 0, so the 0.1 amplitude gate sits well below a single
    copy and admits all gains. Non-focal amplified segments (surrounded by similar-ratio neighbours)
    are filtered as artifacts; focal high-amplitude events are the actionable ones — confirms the
    length/extent notion of "focal."
  - **NCBI Gene — oncogene cytogenetic locations** (rank 5, curated DB) — the registry panel and arms:
    **ERBB2 17q12 · MYC 8q24.21 · EGFR 7p11.2 · CCND1 11q13.3 · MDM2 12q15 · CDK4 12q14.1**. The
    chromosome+arm prefix (17q / 7p / 8q / 11q / 12q) is what a focal amplification's arm label is
    matched against.

- **Documented corner cases / failure modes:** whole-arm event (≥ 98% of arm) excluded as arm-level,
  never focal even when highly amplified; boundary at exactly 0.98 → arm-level (strict `< 0.98` →
  focal); amplitude ≤ `t_amp` (0.1) → not amplified, excluded regardless of length; only focal
  amplifications feed the oncogene mapper; null/empty input handling.

- **Dataset (worked, deterministic — arm length 1,000,000 bp; `t_amp` 0.1; cutoff 0.98):**
  - A (17q, 500 kb, 0.50, log2 1.0) → amplified + focal → **focal amp → ERBB2**.
  - B (8q, 990 kb, 0.99, log2 1.5) → amplified but **arm-level** (0.99 ≥ 0.98) → no.
  - C (7p, 300 kb, 0.30, log2 0.05) → focal but **not amplified** (0.05 ≤ 0.1) → no.
  - D (11q, 980 kb, 0.98, log2 1.0) → **boundary** (0.98 not < 0.98) → arm-level → no.

- **Coverage recommendations (6 items):** MUST — high-amp focal (< 98% arm) reported / ≥ 98% not;
  gain ≤ `t_amp` not reported even if focal; boundary at 0.98 not focal; `IdentifyAmplifiedOncogenes`
  maps 17q→ERBB2, 8q→MYC, 7p→EGFR, 11q→CCND1, 12q→MDM2 **and** CDK4. SHOULD — a non-amplified segment
  is never an oncogene amplification. COULD — null/empty handling.

## Deviations and assumptions

- **ASSUMPTION — amplitude gate combined with the length rule.** GISTIC2.0 classifies focal-vs-broad
  purely by length; the "amplified" gate is the GISTIC2 `t_amp` = 0.1 parameter. Fusing the length rule
  (paper) and the amplitude rule (docs) into one `DetectFocalAmplifications` predicate is this unit's
  **integration choice** — documented, both halves source-backed, not invented.
- **ASSUMPTION — arm fraction supplied as input.** No cytoband table is bundled; the caller supplies
  each segment's arm label + arm length and the algorithm computes `SegLen / ArmLength`. The 0.98
  cutoff and amplitude rule are unchanged.

No source contradictions — Mermel 2011 (98%-of-arm length rule), GISTIC2 docs (`broad_len_cutoff`
0.98 / `t_amp` 0.1), CNVkit (single-copy gain log2 0.585 > 0.1; focal = extent-limited), and NCBI Gene
(oncogene cytobands) corroborate one another across the whole predicate.
