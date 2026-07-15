---
type: source
title: "Evidence: ONCO-FUSION-001 (gene-fusion detection from split + spanning reads, with in-frame codon-phase check)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-FUSION-001-Evidence.md
sources:
  - docs/Evidence/ONCO-FUSION-001-Evidence.md
source_commit: ea13dcc183c950560fe910068244e507f45a455f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-FUSION-001

The validation-evidence artifact for test unit **ONCO-FUSION-001** — **Fusion Gene Detection**
(candidate fusion calling from breakpoint-supporting reads). The **fourteenth ingested unit of
the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized
in [[gene-fusion-detection-read-evidence]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the STAR-Fusion / Arriba split-read + discordant-pair
  + minimum-support paradigm, corroborated across two independent tools and their papers):**
  - **STAR-Fusion source** (Haas et al., rank 3, reference implementation) — the four default
    thresholds: **`MIN_JUNCTION_READS = 1`** (min junction-spanning / split reads),
    **`MIN_SUM_FRAGS = 2`** (min total support = junction reads + spanning frags, requires ≥ 1
    junction read to apply), **`MIN_SPANNING_FRAGS_ONLY = 5`** (min discordant fragments required
    when there are **zero** junction reads), and `MIN_FFPM = 0.1` (abundance filter, recorded but
    not part of the count-threshold rule). Support = **`# junction_reads + # spanning_frags`**.
  - **STAR-Fusion preprint** (Haas et al. 2017 bioRxiv 120295, rank 1; PDF 403, abstract via
    search snippet) — confirms the two evidence classes (discordant + split-read alignments) and
    count-based support scoring "according to the abundance of fusion-supporting reads".
  - **Haas et al. 2019 Genome Biology** benchmark (rank 1) — the 23-method benchmark (incl.
    STAR-Fusion, Arriba) establishes the split-read + discordant-pair + minimum-support paradigm
    as standard.
  - **Arriba output spec** (Uhrig et al., rank 3) — **total support = `split_reads1 + split_reads2
    + discordant_mates`**; a split read's **anchor** is the gene its **longer** aligned segment
    maps to (defines split_reads1 vs split_reads2); **`discordant_mates`** = spanning / bridge
    pairs; the **`reading_frame`** column takes `in-frame` / `out-of-frame` / `stop-codon` / `.`.
  - **Arriba paper** (Uhrig et al. 2021 Genome Research, rank 1) — **split read** = "reads with
    two segments aligning in a noncontiguous fashion"; **discordant mates** (a.k.a. spanning /
    bridge reads) = paired-end reads with mates aligning in a nonlinear way.
  - **Genomics England** gene-fusion reporting (rank 4) + **Wikipedia "Reading frame"** (rank 4,
    for its primary citations Badger & Olsen 1999 / Lodish 6th ed.) — the **exon-phase / modulo-3
    in-frame rule**: a fusion is **in-frame** iff the coding bases the 5' partner contributes
    before the breakpoint, taken **mod 3**, keep the 3' partner's codons in phase; otherwise
    **out-of-frame** and "unlikely to generate a protein".

- **Documented corner cases / failure modes:**
  - **No junction reads** → the `min_sum_frags = 2` rule does not apply; **5** discordant
    fragments required instead (`min_spanning_frags_only`). 1–4 discordant with no split reads →
    filtered (likely artifact).
  - **Junction present but low total** — junction = 1, total = 1 (< 2) → **rejected** (fails
    `min_sum_frags`); junction = 1, total = 2 → **detected** (passes both).
  - **Read-through transcripts** — adjacent same-strand neighbouring genes make common
    false-positive chimeras; guard via support + distinct-gene rules.
  - **Stop codon before junction** — even numerically in-frame, an upstream stop (`stop-codon`
    value) means the 3' partner is not translated (out of scope here — needs transcript
    reconstruction, ONCO-FUSION-003).
  - **Same gene both sides** — `gene5p == gene3p` is not a fusion (registry invariant
    `gene5p ≠ gene3p`) → **rejected**.

- **Datasets (deterministic, derived from the cited threshold rules — no raw BAM table is openly
  downloadable, so candidate-level counts are constructed directly from the rules):**
  - **Breakpoint-evidence candidates** (`split1, split2, discordant`; total = sum; junction =
    split1+split2): **C1 EML4-ALK** (3,2,4 → DETECTED, junc 5 ≥ 1 & sum 9 ≥ 2); **C3 TMPRSS2-ERG**
    (1,0,1 → DETECTED, junc 1 & sum 2); **C4 CD74-ROS1** (0,0,5 → DETECTED, no junc & span 5 ≥ 5);
    **C5 NCOA4-RET** (0,0,4 → REJECTED, span 4 < 5); **C6 KIF5B-RET** (1,0,0 → REJECTED, junc 1 but
    sum 1 < 2); **C7 ALK-ALK** (5,5,5 → REJECTED, gene5p == gene3p). (C2 BCR-ABL1 1,0,0 same as C6
    → REJECTED for sum < 2.)
  - **Reading-frame phase** (`5' coding bases before breakpoint`, `3' start phase`,
    `(5pBases − 3pPhase) mod 3` → frame): **F1** 300, 0 → 0 → **in-frame**; **F2** 301, 0 → 1 →
    out-of-frame; **F3** 302, 0 → 2 → out-of-frame; **F4** 301, 1 → 0 → **in-frame**.

- **Coverage recommendations (9 items):** MUST — junction ≥ 1 AND total ≥ 2 → DETECTED; 0 junction
  → DETECTED only if discordant ≥ 5 (4 rejected); junction = 1 but total = 1 → REJECTED; total =
  split1+split2+discordant; `gene5p == gene3p` → REJECTED; in-frame iff (5' bases − 3' phase) mod 3
  == 0. SHOULD — results ordered by descending total support; null → ArgumentNullException /
  negative counts → ArgumentException. COULD — custom thresholds override the defaults.

## Deviations and assumptions

- **ASSUMPTION — candidate-level input granularity.** The unit consumes already-grouped
  breakpoint candidates with per-class counts (`split_reads1`, `split_reads2`, `discordant_mates`),
  **not** raw BAM records (extracting chimeric reads is a separate `FindChimericReads` method, out
  of the canonical-threshold scope). Mirrors the Arriba output schema; makes the threshold rule
  deterministically testable.
- **ASSUMPTION — in-frame uses coding-base phase, not stop-codon scanning.** In-frame status is
  `(5pCodingBases − 3pStartPhase) mod 3 == 0` per the exon-phase rule; the unit does **not** scan
  the spliced transcript for premature stops (Arriba's `stop-codon` value) — that needs transcript
  reconstruction (ONCO-FUSION-003 scope).

No source contradictions — STAR-Fusion and Arriba agree on the split-read + discordant-pair +
minimum-support model and the additive total-support definition; the in-frame rule is the standard
exon-phase / modulo-3 definition cross-checked against Wikipedia's primary citations.
