---
type: source
title: "Validation report: ANNOT-GENE-001 (prokaryotic gene prediction + Shine-Dalgarno RBS detection)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-GENE-001.md
sources:
  - docs/Validation/reports/ANNOT-GENE-001.md
source_commit: fc49476951d4c26bf663220b231007d8327e59cf
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-GENE-001

The two-stage **validation write-up** for test unit **ANNOT-GENE-001** (ORF-based prokaryotic
gene prediction + Shine-Dalgarno ribosome-binding-site detection on both strands), validated
2026-06-24. This is the *report* artifact that feeds one row of the [[validation-ledger]]; it
records the validator's **verdict** on both the algorithm description and the shipped code. The
algorithm itself is summarized in [[prokaryotic-gene-prediction-rbs]]; the two-stage methodology
is the [[validation-protocol]]. Distinct from any pre-implementation `annot-gene-001-evidence`
artifact.

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** An **independent re-validation** (fresh
context): all consensus/spacing values were re-retrieved from primary sources and every
reverse-strand coordinate was re-derived by hand in Python without lifting the repo's expected
values. **No defect found.** Filtered suite `GenomeAnnotator_Gene_Tests` → **39 passed / 0
failed**. No code or test change required; no follow-up.

## Canonical methods validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`:

- `PredictGenes(dna, minOrfLength, prefix)` — ORF-based gene prediction on both strands; labels
  CDS and sets `protein_length = ProteinSequence.TrimEnd('*').Length` (excludes the stop → `(End−Start)/3 − 1`).
- `FindRibosomeBindingSites(...)` — forward-strand SD only (legacy; `.Where(o => !o.IsReverseComplement)`).
- `FindRibosomeBindingSitesBothStrands(...)` — forward + reverse passes, reports `(position, sequence, score, strand)`.
- private `ScanStrandForShineDalgarno(...)` — the single shared single-strand scan reused by both passes.

## Stage A — description (algorithm faithfulness)

Confirmed against primary sources (re-fetched this session):

- **Shine-Dalgarno consensus `AGGAGG`** (Wikipedia "Shine-Dalgarno sequence"); the shorter code
  motifs `GGAGG / AGGAG / GAGG / AGGA` are all substrings of the consensus, biologically justified
  (T4 early genes are GAGG-dominated).
- **Anti-SD 3' tail of 16S rRNA `5'-YACCUCCUUA-3'`** (Y = pyrimidine); SD located **~8 bases
  upstream** of the AUG.
- **Optimal aligned spacing 5 nt** (Chen et al. 1994, *NAR* 22(23):4953-4957, PMID 7528374), where
  aligned spacing = SD 3' end → the A of the AUG.
- The default `[4,15]` spacing window is an **honestly-declared permissive heuristic search
  window** (the method reports candidates; it does not score optimality) — not a description error.

**Model:** scan all 3 frames on both strands for start→stop, filter by minimum length, label CDS.
The SD is a feature of the mRNA — for a reverse-strand gene it is the reverse complement of the
forward genomic strand, so the reverse-strand SD is found by scanning the **reverse complement**,
not the forward strand.

## Stage B — implementation (the reverse-strand coordinate mapping)

The load-bearing correctness point is the reverse-strand SD position mapping, **independently
re-derived (not lifted from repo tests)**:

- Reverse ORFs are un-mapped back into revComp space (`o with { Start = len−o.End, End = len−o.Start }`),
  scanned on `revComp`, and the motif hit is mapped to forward coordinates via
  **`forwardPosition = len − hit.position − motifLen`** (`GenomeAnnotator.cs:341`).
- Worked derivation: `CreateSequenceWithReverseStrandSd("AGGAGG", spacing 8)` builds
  `S = C×10 + AGGAGG + C×8 + (ATG + A×297 + TAA)`, `|S| = 327`; genomic input = `RC(S)`. Code scans
  `RC(RC(S)) = S`, AGGAGG at revPos 10 → `forwardPos = 327 − 10 − 6 = 311`. Independently verified
  `RC(S)[311..317] = "CCTCCT"` (the anti-SD complement on the forward strand); aligned spacing in
  revComp space = `24 − 10 − 6 = 8` ✓. Other cases re-derived: maxDist-15 → 318; combined-reverse → 638.
- Spacing test in `ScanStrandForShineDalgarno`: `distanceToStart = orf.Start − genomicPos − motif.Length ∈ [minDistance, maxDistance]`; score = `motif.Length / 6.0`. Window semantics include `maxDistance`, exclude above it.

## Test quality + independent mutation check

39 tests (32 legacy + R1-R7 reverse). Assertions are exact positions/sequences/scores/strands or
exact emptiness — no bare "no throw". R1 pre-guards the construct geometry
(`sequence.Substring(311,6) == "CCTCCT"`) so a position assertion cannot silently track a helper
change; R7 asserts byte-for-byte parity of `(position, sequence, score)` between the legacy
forward-only path and the both-strands `+` hits (no forward-behaviour drift). The validator's
**independent mutation** — replacing the reverse mapping with `forwardPosition = hit.position` —
**failed 3 tests (R1/R3/R4)**, confirming the mapping is genuinely locked; source reverted, suite
re-run green 39/39.

## Findings

- **No findings.** Both stages PASS, State ✅ CLEAN. Consensus, anti-SD tail, ~8 nt location, 5 nt
  optimal aligned spacing, the reverse coordinate mapping, boundary/no-site/both-strand and
  protein-length cases all verified against external sources and hand re-derivation. Zero code
  change.
