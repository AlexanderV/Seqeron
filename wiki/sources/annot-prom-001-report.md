---
type: source
title: "Validation report: ANNOT-PROM-001 (prokaryotic promoter motif detection — GenomeAnnotator.FindPromoterMotifs)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-PROM-001.md
sources:
  - docs/Validation/reports/ANNOT-PROM-001.md
source_commit: 7323bb6a053866cc257942f0e337f7990122c62e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-PROM-001

The two-stage **validation write-up** for test unit **ANNOT-PROM-001** (prokaryotic
promoter motif detection — the −10 / −35 box scan), validated 2026-06-24. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the
validator's **verdict** on both the algorithm description and the shipped code. The
algorithm itself is summarized in [[promoter-detection]]; the two-stage methodology is the
[[validation-protocol]] under [[validation-and-testing]]. Distinct from any
pre-implementation `annot-prom-001-evidence` artifact (none exists), and from the
catalog-scan sibling [[regulatory-element-detection]] (MOTIF-REGULATORY-001,
`GenomicAnalyzer.FindMotif`), which reports the same boxes only as *exact hexamer* hits
with no partial-variant scoring.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN — no defect, no code or test change.**
Promoter-motif filter → **20 passed / 0 failed** (28 `[TestCase]`s). No code changed. The
two divergences from the literature (17 bp spacing **not enforced**; **exact-substring** not
mismatch-tolerant/PWM matching) are explicitly declared scope limits of a motif-scan
primitive, not errors.

## Canonical method validated

`GenomeAnnotator.FindPromoterMotifs(string dnaSequence)` →
`IEnumerable<(int position, string type, string sequence, double score)>`
(`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:584-638`).

## Stage A — description (algorithm faithfulness)

Confirmed against **Wikipedia "Promoter (genetics)"** and **"Pribnow box"** (both fetched
live), with **Pribnow (1975)** PNAS 72(3):784-788 and **Harley & Reynolds (1987)** NAR
15(5):2343-2361 as the cited primary refs:

| Property | Validated value |
|---|---|
| −35 box consensus | **TTGACA** (label correct, not swapped) |
| −10 box consensus (Pribnow) | **TATAAT** (AT-rich; `TATAAT` not `TATATA`) |
| Optimal spacer | **17 bp** between boxes (up to ~600-fold strength effect) — noted, **not enforced** |
| −35 per-position probabilities | T .69 / T .79 / G .61 / A .56 / C .54 / A .54 (sum **373**) |
| −10 per-position probabilities | T .77 / A .76 / T .60 / A .61 / A .56 / T .82 (sum **412**) |
| Score formula | `Σ(matched-position p) / Σ(all-6 consensus p)` → ∈ [0, 1] |
| Coordinate convention | authoritative = distance upstream of TSS (negative); impl `position` is a **0-based string index** (substring-scan primitive, not a TSS-anchored predictor) |

The per-position probabilities hard-coded in source match the Wikipedia "Promoter
(genetics)" frequency tables **verbatim**. The report re-derived all **8** score constants
by hand (Python) — the four −35 variants and four −10 variants — and every value matched the
spec table and code constants exactly (e.g. `TTGAC` = (69+79+61+56+54)/373 = **0.855**;
`ATAAT` = (76+60+61+56+82)/412 = **0.813**; full boxes = **1.000**). Both divergences (no
spacing enforcement; exact-substring rather than mismatch-tolerant matching) are declared in
the TestSpec "Deviations from Literature" section — legitimate scoping choices for a
motif-scan primitive, so **Stage A: PASS**.

## Stage B — implementation

Motif tables `minus35Motifs` / `minus10Motifs` hold the full 6-mer plus **prefix-5,
suffix-5, prefix-4** variants with the literature-derived scores above — reproduced directly
from the code's constant arrays and all 8 agree with the Stage-A hand computation. Input is
upper-cased via `ToUpperInvariant()` → case-insensitive (M07). Each motif scanned with
`for (i = 0; i <= seq.Length - motif.Length; i++)` and an exact `Substring` compare; every
occurrence yielded with its 0-based `position`. Loop bound is safe for empty / short input
(`motif.Length > len` ⇒ no iterations) so empty and all-C sequences yield nothing (M05/M06).
Labels correct: `-35 box` ↔ TTGACA, `-10 box` ↔ TATAAT (not swapped). Single canonical static
method — no `*Fast` variant or instance delegate to cross-check. Edge cases traced:
no-motif → empty; one box only → only that box's variants; multiple same-type → all positions
(M08: TTGACA at 0 and 11). **Stage B: PASS.**

## Test quality

`GenomeAnnotator_PromoterMotif_Tests.cs`: 20 methods → 28 `[TestCase]`s. Assertions check
**exact** sourced values — positions, `type` strings, scores to ±0.001 — plus monotonic
ordering (S02) and exact variant decomposition (S01: `TATAAT` → variants at 0,0,1,0). Not
tautological. Real-world case C01 verifies both boxes with ~18 bp spacing at exact positions
4 and 28.

## Findings

- **No defect.** Both stages pass; End state ✅ CLEAN, zero code/test change.
- **Declared scope limits (not defects):** 17 bp `-35`/`-10` spacing is **not** validated
  (independent scans, no promoter-pair assembly); matching is **exact** against a fixed
  8-variant library (no mismatch tolerance, no PWM/HMM). Downstream filtering must supply
  spacing/pairing.
- **`position` is a 0-based string index**, not the biological TSS-relative negative
  coordinate — a substring-scan primitive, consistently scoped and documented.
- **No null guard:** the method dereferences `dnaSequence` immediately in `ToUpperInvariant()`;
  empty string is safe (no valid start positions) but null is not explicitly handled.
- **Score-table assumption (accepted):** the −10 constants follow the repository's "Promoter
  (genetics)" table rather than the alternate Harley-et-al. summary on the "Pribnow box" page;
  ANNOT-PROM-001 standardizes on the chosen table (internally consistent across code + tests).
