---
type: source
title: "Validation report: CHROM-KARYO-001 (karyotype analysis — AnalyzeKaryotype + DetectPloidy)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-KARYO-001.md
sources:
  - docs/Validation/reports/CHROM-KARYO-001.md
source_commit: fcb5a4bc5a9d45258189b05e8014c648d0c5197e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-KARYO-001

The two-stage **validation write-up** for test unit **CHROM-KARYO-001** (karyotype analysis —
karyotyping a set of chromosome descriptors with `AnalyzeKaryotype`, and estimating whole-genome ploidy
from read depth with `DetectPloidy`), validated 2026-06-24. This is the *report* artifact that feeds one
row of the [[validation-ledger]] under the [[validation-protocol]] two-stage process; the wider campaign
is [[validation-and-testing]]. The two algorithms, their invariants, the documented cytogenetic oracles
and the five design decisions are synthesized in the concept [[karyotype-analysis]] (the chromosome-count
karyotyping + ploidy anchor); [[test-unit-registry]] defines the unit. Distinct from
[[chrom-karyo-001-evidence]] — the pre-implementation evidence artifact sourced from `docs/Evidence/` —
this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No code changed; no test change required. All 36
unit tests in `ChromosomeAnalyzer_Karyotype_Tests.cs` assert exact values and pass.

## Scope clarification (important — resolves the prompt's framing)

The validation session was framed around the **Levan, Fredga & Sandberg (1964) arm-ratio /
centromeric-index** chromosome classification (metacentric / submetacentric / subtelocentric /
telocentric). The report explicitly finds that model is **NOT CHROM-KARYO-001**: per
`ALGORITHMS_CHECKLIST_V2.md` lines 73–74 the Levan classification is a separate unit, **CHROM-CENT-001**
(`CalculateArmRatio` / `ClassifyChromosomeByArmRatio` / `AnalyzeCentromere`, synthesized in
[[centromere-analysis]] / [[chrom-cent-001-report]]). CHROM-KARYO-001 is the **karyotype / ploidy /
aneuploidy** unit: `AnalyzeKaryotype` + `DetectPloidy`, sourced to Wikipedia Karyotype / Ploidy /
Aneuploidy. This report validates the actual karyotype unit as defined by its TestSpec, Evidence doc,
checklist block, and test file — the arm-ratio framing is a follow-up on CHROM-CENT-001.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`:

- `AnalyzeKaryotype(chromosomes, expectedPloidyLevel)` (`:136–180`) — separates sex chromosomes from
  autosomes, groups autosomes by base name via `GetChromosomeBaseName` (`:200–212`, strips a trailing
  numeric `_N` copy suffix only), counts copies per group, and flags a group abnormal iff
  `count != expectedPloidyLevel`. Summary fields: `TotalGenomeSize = Σ lengths`,
  `MeanChromosomeLength = total/count`, sex/autosome split by `IsSexChromosome`.
- `GetAneuploidyTerm(count)` (`:186–195`) — static label map 0→Nullisomy, 1→Monosomy, 2→Disomy,
  3→Trisomy, 4→Tetrasomy, 5→Pentasomy, ≥6→"Polysomy (N copies)", matching the Wikipedia Aneuploidy
  table verbatim on an **absolute copy-count** basis.
- `DetectPloidy(normalizedDepths, expectedDiploidDepth)` (`:217–241`) — `ratio = median(depths) /
  expectedDiploidDepth`; `ploidy = round(ratio×2)` clamped to `[1, 8]`;
  `confidence = max(0, 1 − |ratio×2 − ploidy|×2)`. Uses a **true median** (averages the two middle
  elements for even counts, `:227–229`).
- Supporting `Karyotype` record (`:45–53`).
- Tests: `tests/.../ChromosomeAnalyzer_Karyotype_Tests.cs` (36 tests).

## Stage A — description (algorithm faithfulness)

Confirmed against two Wikipedia sources fetched 2026-06-24:

- **Aneuploidy** — the copy-count ladder (Nullisomy 0 / Monosomy 1 / Disomy 2 / Trisomy 3 / Tetrasomy 4
  / Pentasomy 5) is confirmed as **absolute, not relative to expected ploidy**: disomy is normal for
  diploids but aneuploid "in organisms that normally have three or more copies". So in a tetraploid
  context (expected = 4) the code correctly labels 2 copies → "Disomy" and 3 copies → "Trisomy".
- **Ploidy** — haploid 1 / diploid 2 / triploid 3 / tetraploid 4 sets; human 2n = 46, n = 23.

`DetectPloidy` is a **depth-ratio heuristic** (not a literature formula but a reasonable engineering
model): diploid ratio ≈ 1 → ploidy 2, tetraploid ratio ≈ 2 → ploidy 4, agreeing with the ploidy
definition. Summary quantities (total genome size, mean chromosome length, sex/autosome split) are
standard karyotype fields.

Edge-case semantics documented and tested: empty chromosomes → empty karyotype, no aneuploidy (DD1);
empty depths → `(2, 0)` default diploid, zero confidence (DD2); ploidy clamp `[1, 8]` (DD3); Nullisomy
(0 copies) unreachable via `GroupBy` — absent chromosomes form no group, term mapped for completeness
(DD4).

**Independent hand cross-check of `DetectPloidy` (7 rows, all reproduced against test assertions):**

| ratio | round(ratio×2) | clamp | confidence = max(0, 1 − \|ratio×2 − ploidy\|×2) |
|-------|----------------|-------|--------------------------------------------------|
| 1.0 | 2 | 2 | 1.0 |
| 2.0 | 4 | 4 | 1.0 |
| 0.5 | 1 | 1 | 1.0 |
| 1.5 | 3 | 3 | 1.0 |
| 1.2 | 2 | 2 | 1 − 0.4×2 = 0.2 |
| 10.0 | 20 | 8 | max(0, 1 − 12×2) = 0 |
| 0.1 | 0 | 1 | max(0, 1 − 0.8×2) = 0 |

**Finding (internal-consistency note, not a defect):** confidence uses the **clamped** ploidy as the
reference, so for out-of-range depths the |·| term is the distance to the clamp (→ 0 confidence) rather
than to the raw round — internally consistent and explicitly tested (`ExtremeHighDepth` /
`ExtremeLowDepth`). **Findings: none affecting the unit.**

## Stage B — implementation

The aneuploidy-term switch matches the sourced table verbatim; the abnormality trigger is
`count != expectedPloidyLevel`; `DetectPloidy` realises the documented ratio / round / clamp / confidence
steps with a true median. The 7 hand computations above reproduce the in-test comments exactly, and the
36 tests assert exact values (e.g. `Abnormalities[0] == "Disomy chr1"`, confidence `0.2 .Within(1e-10)`,
ploidy/confidence pairs) — covering all M/S/C cases, both invariant sets
(`Total = Autosomes + Sex`; `HasAneuploidy ↔ non-empty Abnormalities`; `Ploidy ∈ [1,8]`;
`Confidence ∈ [0,1]`), and every edge case (empty inputs, single value, clamps, between-ploidy,
absolute-terminology in a tetraploid context). Real, deterministic, value-exact — no tautological
no-throw assertions.

## Boundary vs sibling units

- **vs CHROM-ANEU-001 / [[aneuploidy-detection]]:** a sibling method `IdentifyWholeChromosomeAneuploidy`
  (`:881–917`) reuses the same term vocabulary but **hardcodes a diploid baseline** (`CopyNumber != 2`,
  no `2 => "Disomy"` case). That is correct *for that method's diploid-only contract* and belongs to
  CHROM-ANEU-001 (`DetectAneuploidy`, [[chrom-aneu-001-report]]) — not a divergence within
  CHROM-KARYO-001. KARYO owns the descriptor-count karyotype + depth-ratio ploidy path; ANEU owns the
  read-depth→CN path. No overlap conflict.
- **vs CHROM-CENT-001 / [[centromere-analysis]]:** the Levan arm-ratio classification (see scope
  clarification above) is CENT, not KARYO.

## Findings & follow-ups

- **No code defect and no test change.** Terminology (absolute copy count), the ploidy/depth-ratio
  mapping, and the seven `DetectPloidy` rows were independently confirmed against Wikipedia
  "Aneuploidy" / "Ploidy" and reproduced by hand; the code computes the validated formula exactly and
  the 36 tests lock the sourced values.
- **Process follow-up:** the prompt's Levan arm-ratio framing applies to **CHROM-CENT-001**; validate
  that separately (already covered by [[chrom-cent-001-report]]).
