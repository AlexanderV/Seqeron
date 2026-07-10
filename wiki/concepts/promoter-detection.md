---
type: concept
title: "Promoter detection (scored −10 / −35 box motif scan)"
tags: [motif, algorithm, annotation]
sources:
  - docs/Validation/reports/ANNOT-PROM-001.md
  - docs/algorithms/Annotation/Promoter_Detection.md
source_commit: 7323bb6a053866cc257942f0e337f7990122c62e
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: alternative_to
      object: concept:regulatory-element-detection
      source: annot-prom-001-report
      evidence: "Both detect the prokaryotic −10 (TATAAT) and −35 (TTGACA) boxes but differ in method and output: promoter-detection (GenomeAnnotator.FindPromoterMotifs) reports scored hits with partial-variant decomposition (full 6-mer + prefix-5/suffix-5/prefix-4, score = Σ matched-position p / Σ all-6 p ∈ [0,1]); regulatory-element-detection (GenomicAnalyzer.FindMotif catalog scan, MOTIF-REGULATORY-001) reports the same boxes only as exact hexamer hits with no scoring."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-prom-001-report
      evidence: "Test Unit ID: ANNOT-PROM-001 (docs/algorithms/Annotation/Promoter_Detection.md), validated 2026-06-24, Stage A/B both PASS, End state CLEAN."
      confidence: high
      status: current
---

# Promoter detection (scored −10 / −35 box motif scan)

**Promoter detection** here is `GenomeAnnotator.FindPromoterMotifs(string dnaSequence)`
(`GenomeAnnotator.cs:584-638`): a **bacterial motif-search helper**, not a full
transcription-start-site predictor. It scans DNA for the prokaryotic **−35 box (`TTGACA`)**
and **−10 box / Pribnow box (`TATAAT`)** consensus substrings plus a small set of derived
prefix/suffix variants, and assigns each match a **literature-based score** derived from
E. coli position-specific nucleotide occurrence frequencies. It returns
`(int position, string type, string sequence, double score)` per hit. Validated under test
unit **ANNOT-PROM-001**; the validation write-up is [[annot-prom-001-report]], and the
two-stage methodology sits under [[validation-protocol]] / [[validation-and-testing]].

## Scoring model

For a matched consecutive substring of one consensus element the score is the fraction of
the full-consensus probability weight that the matched positions carry:

`score = Σ(pᵢ for matched consensus positions) / Σ(pᵢ for all 6 consensus positions)`

with the per-position E. coli frequencies (matching the Wikipedia "Promoter (genetics)"
table verbatim):

- **−35 `TTGACA`:** T .69 / T .79 / G .61 / A .56 / C .54 / A .54 → denominator **373**
- **−10 `TATAAT`:** T .77 / A .76 / T .60 / A .61 / A .56 / T .82 → denominator **412**

Invariants: every score ∈ **[0, 1]** (INV-01); a full `TTGACA` / `TATAAT` match scores
**1.0** (INV-02); every reported motif belongs to the fixed variant library (INV-03).

## The fixed 8-variant library

A single full-consensus occurrence emits **multiple overlapping hits** — the full 6-mer plus
prefix-5, suffix-5, and prefix-4 variants are all reported separately when present. All eight
score constants were independently hand-recomputed (Python) during validation and match code
and spec exactly:

| Class | Variant | Score | | Class | Variant | Score |
|---|---|---|---|---|---|---|
| −35 box | `TTGACA` | 1.000 | | −10 box | `TATAAT` | 1.000 |
| −35 box | `TTGAC` | 0.855 | | −10 box | `TATAA` | 0.801 |
| −35 box | `TGACA` | 0.815 | | −10 box | `ATAAT` | 0.813 |
| −35 box | `TTGA` | 0.710 | | −10 box | `TATA` | 0.665 |

## Contract and edge cases

| Aspect | Behaviour |
|---|---|
| `position` | **0-based index into the input string** — NOT the biological TSS-relative negative coordinate. A substring-scan primitive, not a TSS-anchored predictor. |
| `type` | `-35 box` or `-10 box` (labels correct, not swapped) |
| Case | input `ToUpperInvariant()`-ed first → lowercase/mixed-case accepted (M07) |
| Empty sequence | no hits (scan loops never execute) |
| Null sequence | **not guarded** — dereferenced immediately in `ToUpperInvariant()` |
| Multiplicity / overlap | all supported hits reported; −35 and −10 scanned **independently** |
| Complexity | `O(n × v)`, `v = 8` variants, motif lengths 4–6; `O(1)` + output |

## Declared scope limits (not defects)

Two divergences from the literature are **deliberate** and locked by the TestSpec, not hidden
defects:

- **No `-35`/`-10` spacing check.** The canonical ~**17 bp** optimal spacer (up to ~600-fold
  strength effect) is noted in the sources but **not enforced** — the method never assembles
  complete promoter pairs. Downstream filtering must supply pairing/spacing.
- **Exact-substring matching only.** Real promoters tolerate 2–3 mismatches; this library
  matches exact substrings with no mismatch tolerance and no PWM/HMM/sigma-factor model. It is
  therefore **motif annotation**, not full promoter prediction.

The −10 score table follows the repository's "Promoter (genetics)" values rather than the
alternate Harley-et-al. summary on the "Pribnow box" page — an accepted assumption that
ANNOT-PROM-001 standardizes on (internally consistent across code and tests).

## Scope and siblings

This is the **scored −10/−35 detector**. It is the alternative to the
catalog-scan [[regulatory-element-detection]] (MOTIF-REGULATORY-001,
`GenomicAnalyzer.FindMotif`), which reports the same two boxes — alongside TATA/CAAT/GC/Kozak/
Shine-Dalgarno/poly(A)/TF sites — only as **exact hexamer** hits with **no partial-variant
scoring**. Where this detector emits scored 4–6-base variants and their score fractions, the
catalog scanner emits bare `Name`/`Pattern`/`Sequence` occurrences. Both are exact-match
primitives (no PWM), distinct from the position-weight-matrix scoring in the splicing family
(cf. [[splice-acceptor-site-prediction]]) and from the strand/spacing-aware ribosome-binding-
site finder [[prokaryotic-gene-prediction-rbs]]. **No source contradictions.**
