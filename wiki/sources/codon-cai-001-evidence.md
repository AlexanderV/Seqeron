---
type: source
title: "Evidence: CODON-CAI-001 (Codon Adaptation Index — CAI)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-CAI-001-Evidence.md
sources:
  - docs/Evidence/CODON-CAI-001-Evidence.md
source_commit: 06dbe159ef5ed5e03f98988579713403bc50f51b
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-CAI-001

The validation-evidence artifact for test unit **CODON-CAI-001** (Codon Adaptation Index — CAI
calculation). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is summarized in
[[codon-adaptation-index]], which builds on the [[relative-synonymous-codon-usage]] family
normalization. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** — Wikipedia "Codon Adaptation Index" for the definition (most widespread
  codon-bias technique; deviation of a gene from a reference gene set; expression-level predictor)
  and the `w_i = f_i/max(f_j)` + `CAI = (∏ w_i)^{1/L}` formulae; **Sharp & Li (1987)**, *Nucleic
  Acids Res.* 15(3):1281–1295 (PMID 3547335) as the original CAI paper; **Jansen, Bauer & Stadler
  (2003)** (PMC2684136), quoted **verbatim** for the single-codon-amino-acid exclusion rule and its
  rationale (Met/Trp `w≡1` regardless of bias inflates CAI for Met/Trp-rich genes); Kazusa Codon
  Usage Database as the reference-table source.
- **Algorithm spec** — geometric mean of relative adaptiveness `w_i` (codon frequency ÷ most-used
  synonymous codon's frequency); equivalently `exp((1/L)·Σ ln w_i)`. Range [0,1]; CAI=1 all-optimal,
  →0 as rare codons accumulate; stop codons excluded; single-codon Met/Trp always `w=1`.
- **Dual convention** — default `CalculateCAI(seq, table)` **includes** Met/Trp (`w=1`, historical);
  `CalculateCAI(seq, table, excludeSingleCodonAminoAcids: true)` **excludes** them per the canonical
  Sharp & Li 1987 / Jansen 2003 rule (excluding can drop L to 0 → CAI 0).
- **Datasets** — E. coli K12 (Kazusa 316407) Leu and Arg codon tables with hand-computed `w`;
  worked cases `AUG`→1.0, `CUG-CCG-ACC`→1.0, `CUA-CCA-ACA`→0.1980; exclusion-mode cases
  `AUGUGG`→0 (both excluded, L=0), `AUGCUACUA` incl 0.18566…/excl 0.08, `AUGUGGCUA` incl 0.43088…/
  excl 0.08, `CUGCUA` 0.28284… (flag no-op).
- **Edge cases** — empty sequence → 0; Met/Trp-only → 1.0 (inclusive); all-optimal → 1.0; all-rare
  → ~0; same sequence differs by organism (codon usage varies).

## Assumptions and the one deviation (from the artifact)

1. **`1e-6` clamp deviation.** When a codon's frequency is 0 but its amino acid has other codons in
   the table (`maxFreq > 0`), `w` is clamped to `1e-6` instead of 0 — an incomplete-codon-table
   protection so one missing entry cannot zero the whole gene's CAI (strict Sharp & Li would give
   `w=0` → `log(0)`). Unknown amino acid or `maxFreq=0` → `w=NaN`, codon skipped by caller.
2. **Empty sequence → 0** by convention (no codons to evaluate).
3. Frequency tables verified against the Kazusa database (March 2026); implementation converts
   T→U, is case-insensitive, and splits into codons before scoring.

**Contradictions:** none between sources — Wikipedia's formulae, Sharp & Li 1987, and the Jansen
2003 verbatim exclusion quote agree. **One cross-page nuance to reconcile** (not a source conflict):
[[relative-synonymous-codon-usage]] describes CAI's zero-codon guard as a **0.5 pseudocount**
(Sharp & Li's reference-table convention), whereas this artifact documents the Seqeron
**implementation** using a **`1e-6` clamp** at score time — different values, different stage, both
avoiding `log(0)`. Flagged on [[codon-adaptation-index]].
