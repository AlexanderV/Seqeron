---
type: source
title: "Validation report: CHROM-SYNT-001 (synteny analysis — collinear blocks + rearrangement detection, ChromosomeAnalyzer.FindSyntenyBlocks / DetectRearrangements)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-SYNT-001.md
sources:
  - docs/Validation/reports/CHROM-SYNT-001.md
source_commit: 7a7cdd292084b683c3ee8baa98bbbcd61e441c4b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-SYNT-001

The two-stage **validation write-up** for test unit **CHROM-SYNT-001** (synteny analysis — building
collinear synteny blocks from ortholog pairs and classifying the chromosomal rearrangements between
adjacent blocks), validated 2026-06-24. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The two algorithms, their invariants, the documented oracles and the
edge cases are synthesized in the concept [[synteny-and-rearrangement-detection]] (the shared synteny
anchor for the chromosome- and comparative-genomics families); [[test-unit-registry]] defines the
unit. Distinct from [[chrom-synt-001-evidence]] — the pre-implementation evidence artifact sourced
from `docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN.** No code defect; no code change
required. The synteny filter (`FullyQualifiedName~ChromosomeAnalyzer_Synteny_Tests`) ran **19 passed,
0 failed**; every asserted count, coordinate, strand and size reproduced. This is an independent
re-validation from a fresh context — the source has **not changed** since the prior archived report
(`git diff cb113ce..HEAD` over the synteny region is empty; the last code touch `f0568862` only
strengthened tests), so conclusions were re-derived rather than copied.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (synteny region
lines 638–823):

- `FindSyntenyBlocks(orthologPairs, minGenes=3, maxGap=10)` (`:643–716`) — groups pairs by
  `(Chr1, Chr2)`; sorts by `Start1`; walks consecutive pairs with `currentForward = curr.Start2 >
  prev.End2`. The first step (i==1) fixes block direction `isForward`; a pair stays collinear iff
  `currentForward == isForward` **and** `gap1 ≤ maxGap·1e6` **and** `gap2 ≤ maxGap·1e6`
  (`gap2 = Math.Abs(...)`). On break, or at the last index, emits a block when `geneCount ≥ minGenes`;
  strand `= isForward ? '+' : '-'`; Species2 coords use `Min(first,last)` / `Max(first,last)` to stay
  valid under inversion; `SequenceIdentity = NaN`.
- `DetectRearrangements(syntenyBlocks)` (`:721–823`) — sorts blocks by `(Chr1, Start1)`; for adjacent
  same-`Chr1` pairs — different `Chr2` → **Translocation**; same `Chr2` + different strand →
  **Inversion**; same `Chr2` + same strand + asymmetric gap (`gap1 > 0 && gap2 ≥ 0 && gap1 > 2·gap2`)
  → **Deletion** (Size `= gap1 − max(0,gap2)`). A separate O(n²) pass flags overlapping Species1
  regions mapping to different Species2 targets → **Duplication**.
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Synteny_Tests.cs` (19 tests).

## Stage A — description (algorithm faithfulness)

Confirmed against **MCScanX** (Wang et al. 2012, *Nucleic Acids Res.* 40(7):e49, re-read via
PMC3326336), **Wikipedia "Synteny"** and **Wikipedia "Chromosomal rearrangement"**. A collinear
(synteny) block is a chain of pairwise-collinear anchor genes preserving gene order between two
genomes; MCScanX's two parameters **MATCH_SIZE** (min anchors, default 5) and **MAX_GAPS** (default
20) map directly onto the spec's `minGenes` and `maxGap`. Forward (`+`) blocks have target positions
increasing with reference; inverted (`-`) blocks have them decreasing; a gap over threshold or an
order reversal breaks the run; runs with `< minGenes` anchors are dropped. `SequenceIdentity = NaN`
is sound because identity is not computable from coordinate-only input (MCScanX derives it from BLAST
alignments). The four rearrangement types — Inversion, Translocation, Deletion, Duplication — match
the encyclopedic definitions.

**Independent hand cross-check** (reference order 1,2,3,4): target `1,2,3,4` → 1 forward block of 4
(`+`); `4,3,2,1` → 1 inverted block of 4 (`-`); `1,2,5,3,4` → 1 forward block of 3 `{1,2,5}` with the
trailing `{3,4}` dropped (`< minGenes=3`, since the 5→3 step reverses direction). All match spec
M1/M2 and the gap/break semantics.

**Stage A notes (documentation-only, non-blocking):**

1. Strand / inverted-block classification is **not** in the Wikipedia "Synteny" article (which is
   strand-silent), the cited label for oracle M2. The `+`/`-` distinction is nonetheless standard and
   authoritative via **MCScanX dot-plots and SyRI inversion calls** — the definition is correct, only
   the *citation label* is weak.
2. The spec's collinearity rule is the simpler "consecutive consistent ordering," not MCScanX's full
   scoring DP / LIS chaining. Both are valid block definitions; the spec's is a legitimate documented
   simplification.

## Stage B — implementation

Every worked example was re-traced against the code and reproduced exactly: M1 forward (4 genes) → 1
block `+`, Species1/2 1000–8000; M2 reverse (4 genes) → 1 block `-`, Species2 `Min(8000,2000)=2000` /
`Max(9000,3000)=9000`; M5 exactly `minGenes=3` → 1 block `+`, 1000–6000; M6 two chr pairs → 2 blocks
of 3; M16 gap split (`maxGap=2` ⇒ 2 MB; gene3→gene4 gap 2,994,000 bp > 2 MB) → 2 blocks of 3; M9
Inversion → Pos1 50000 / Pos2 60000 / Size 10000; M10 Translocation → Chr2 "chrB", Pos1 50000, Pos2
1000, Size null; M14 Deletion → gap1 100000, gap2 5000, `100000 > 2·5000` → Size 95000; M15
Duplication → overlap 20000–50000, different Species2 → Size 30000. Edge cases (empty→empty; single
marker→empty; single block→empty; fully collinear genome→empty) all hold.

**Numerical robustness:** `maxGap·1e6` is evaluated as `int` and fits in `int.MaxValue` (~2.1e9) for
`maxGap` up to ~2147 (default 10 safe); `gap2` uses `Math.Abs` to guard inverted runs; no div-by-zero.
There are **no `*Fast`/delegate variants** — two independent canonical methods, consistent.

**Test-quality audit:** 19 tests assert exact hand-calculated values (counts, coordinates, strand,
sizes), deterministic, covering forward / inverted / gap-split / min-genes / multi-chromosome and all
four rearrangement types plus invariants (strand ∈ {+,−}, Start ≤ End, NaN identity, valid type
strings, Position1 set). Strong, not tautological.

**Stage B notes (documented simplifications, not defects):** (1) greedy pairwise chaining rather than
DP/LIS — locally-monotonic but globally-interleaved markers are absorbed into a block, faithful to the
spec's definition and unasserted-against; (2) `DetectRearrangements` is adjacency-only with heuristic
thresholds (Deletion uses `gap1 > 2·gap2`) — reasonable coordinate-only heuristics, internally
consistent and tested, not the full SyRI classifier.

## Findings

- **No code defect and no test change (State CLEAN).** Every worked example and cross-check value
  reproduced exactly; divergences are documented simplifications, not defects.
- **Documentation-only follow-ups (non-blocking):** re-attribute the inverted-block (M2) and gap-split
  (M16) tests to **MCScanX / SyRI** rather than the strand-silent Wikipedia "Synteny" page; reconcile
  the spec / Evidence "20 tests" with the actual **19** in the test file (a documentation off-by-one,
  not a coverage gap).
