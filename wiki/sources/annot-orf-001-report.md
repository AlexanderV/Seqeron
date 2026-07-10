---
type: source
title: "Validation report: ANNOT-ORF-001 (annotation-layer ORF detection — GenomeAnnotator.FindOrfs)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-ORF-001.md
sources:
  - docs/Validation/reports/ANNOT-ORF-001.md
source_commit: 23decb5c5e895bf2baa626a971bfb7de3b02322b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-ORF-001

The two-stage **validation write-up** for test unit **ANNOT-ORF-001** (annotation-layer
open-reading-frame detection), validated 2026-06-24. This is the *report* artifact that feeds
one row of the [[validation-ledger]]; it records the validator's **verdict** on both the
algorithm description and the shipped code. The algorithm itself is summarized in
[[open-reading-frame-detection]]; the two-stage methodology is the [[validation-protocol]] under
[[validation-and-testing]]. Distinct from any pre-implementation `annot-orf-001-evidence`
artifact, and from the sibling ATG-only `GenomicAnalyzer` unit [[genomic-orf-001-evidence]].

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End state: ✅ CLEAN — no defect, no code or test
change.** ORF filter → **35 passed / 0 failed**; full project **18208 passed / 0 failed**. The
only note is the **non-canonical `requireStartCodon=false` run-off seeding** path (outside
standard-ORF scope); the canonical `requireStartCodon=true` path is fully correct.

## Canonical methods validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`:

- `FindOrfs(dna, minLength, searchBothStrands, requireStartCodon)` — forward scan then
  reverse-complement scan with coordinate adjustment (`FindOrfsInStrand` frames 0/1/2 →
  `FindOrfsInFrame` start/stop scanning + nested-start handling + run-off handling).
- `FindLongestOrfsPerFrame(dna, searchBothStrands)` — groups by signed frame, picks longest
  protein; reuses `FindOrfs(minLength:1)`, returns 6 signed-frame keys (or 3 forward-only).

## Stage A — description (algorithm faithfulness)

Confirmed against **Rosalind "ORF"** (fetched live), **Wikipedia "Open reading frame"**,
**NCBI ORFfinder**, **Deonier et al. (2005)** (~150 codons) and **Claverie (1997)** (~100 codons):

| Property | Validated value |
|---|---|
| Start codons | ATG canonical **+ GTG, TTG** (alt starts; NCBI option, sourced superset of Rosalind's ATG-only) |
| Stop codons | TAA, TAG, TGA |
| Reading frames | 6 = 3 forward + 3 reverse-complement |
| ORF end | first in-frame stop, no internal in-frame stop |
| Stop in coordinates | **included** in `Sequence`/`End`; `ProteinSequence` ends with `*`, excluded from peptide |
| Coordinate base | 0-based, half-open `[Start, End)` |
| `minLength` unit | **amino acids** (codons before the stop) |
| Nested/overlapping | every start codon upstream of a shared stop yields its own ORF |

Independent cross-check reproduced the authoritative Rosalind sample's **four** distinct proteins
`{MLLGSFRLIPKETLIQVAGSSPCNLS (reverse strand), M, MGMTPRLGLESLLE, MTPRLGLESLLE}`. The report flags
that the orchestration prompt quoted a **truncated** Rosalind sequence; validation used the
authoritative full sequence.

## Stage B — implementation

Live-run cross-verification against the Rosalind full sample matched all four canonical proteins;
the only extras (`LD`, `LE`) arise solely from the sourced GTG/TTG alternative starts, and M12
asserts **containment** so the superset causes no false failure. Load-bearing points verified in
code: start/stop sets `{ATG,GTG,TTG}`/`{TAA,TAG,TGA}` (`OrdinalIgnoreCase`); six frames on `dna` +
`GetReverseComplementString(dna)`; on a stop, all open starts complete and `currentOrfStarts`
clears (no ORF crosses an in-frame stop); reverse-strand mapping `adjStart = len − End`,
`adjEnd = len − Start` maps half-open intervals with no off-by-one; nested ORFs share an `End`
(Rosalind `MGMTPRLGLESLLE` start 24 nests `MTPRLGLESLLE` start 30, shared end 69). Edge cases
(empty/null, no start, no stop, below/exactly `minLength`, boundary ORF ending at seq end, alt
starts, N-in-codon, zero-length run-off guard) all traced to code lines and tests. 35 tests, exact
assertions (counts, positions, frame numbers, protein strings, reverse-strand protein).

## Findings

- **No defect.** Both stages effectively pass; End state ✅ CLEAN, zero code/test change.
- **PASS-WITH-NOTES (minor):** the `requireStartCodon=false` run-off logic does not emit a
  coding stretch from the frame origin before the first start/stop — non-default path, outside
  validated standard-ORF semantics (Rosalind/NCBI require a start codon); the zero-length run-off
  case is correctly guarded (`endPos > start`).
- **Stop-inclusion convention:** `Sequence`/`End` include the stop while the `minLength`
  amino-acid count excludes it; both internally consistent and locked by tests (NCBI-style
  coordinates include the stop).
