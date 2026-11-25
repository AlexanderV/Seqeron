---
type: concept
title: "Open reading frame detection (six-frame ATG→stop enumeration)"
tags: [annotation, algorithm]
mcp_tools:
  - find_open_reading_frames
sources:
  - docs/Evidence/GENOMIC-ORF-001-Evidence.md
  - docs/algorithms/Analysis/Open_Reading_Frame_Detection.md
  - docs/algorithms/Annotation/ORF_Detection.md
  - docs/Validation/reports/ANNOT-ORF-001.md
  - docs/Validation/review-2026-09/B11.md
source_commit: 8b85d3b5496273e86f3be2eef7f720a6ce7a257b
created: 2026-07-09
updated: 2026-11-25
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-orf-001-evidence
      evidence: "Test Unit ID: GENOMIC-ORF-001 ... Algorithm: Open Reading Frame (ORF) Detection — GenomicAnalyzer.FindOpenReadingFrames"
      confidence: high
      status: current
---

# Open reading frame detection (six-frame ATG→stop enumeration)

> **Updated 2026-11-25 (B11 F1–F3, F41).** `GenomeAnnotator.FindOrfs` no longer has its own ORF engine: it delegates to `Translator.FindOrfs` (stop-to-stop, = EMBOSS `getorf -find 0`) / `Translator.FindStartToStopOrfs`; alternative initiators translate as Met, `requireStartCodon=false` keeps leading/stop-less regions, duplicate sub-ORFs are gone.


**ORF detection** enumerates candidate protein-coding spans in a DNA sequence: each
begins at a start codon (**ATG**) and ends at the **first in-frame stop** codon
(TAA/TAG/TGA) with no internal in-frame stop. Because double-stranded DNA can be read in
three frames on each strand, the search covers **six frames** (three forward, three on
the reverse complement). Seqeron exposes it as `GenomicAnalyzer.FindOpenReadingFrames`,
an exact, deterministic linear-time codon scan (no heuristics, no coding-potential
model). Validated under test unit **GENOMIC-ORF-001**; the validation record is
[[genomic-orf-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Core model

For strand `S` and each frame offset `f ∈ {0,1,2}`, codons are `S[f..f+3), S[f+3..f+6), …`.
An ORF is a span `S[a..b+3)` where `S[a..a+3) = ATG`, `S[b..b+3)` is a stop, `(b−a) mod 3 = 0`,
and no codon strictly between `a` and `b` is a stop. The **reported span includes the stop
codon** (so `Length % 3 == 0`); the **translated protein candidate excludes** it
(translation runs "until a stop", Rosalind). The reverse complement is scanned identically,
marking results `IsReverseComplement = true`. Standard [[genetic-code-translation|genetic code]],
NCBI transl_table=1 (start ATG; stops TAA/TAG/TGA).

## Nested ORFs sharing a stop (the correctness rule)

Every in-frame ATG that reaches a downstream in-frame stop is reported — including a
**downstream ATG in the same frame that shares the same stop**, which opens a second,
shorter ORF. Both protein candidates are returned. This is the canonical Rosalind
semantics (sample yields both `MGMTPRLGLESLLE` and the nested `MTPRLGLESLLE`). A greedy
scan that stops at the first ATG per stop under-reports and is a defect (a pre-existing
greedy bug was fixed to this behavior). An ATG with **no** downstream in-frame stop is an
incomplete ORF and yields **no** candidate.

## API contract and invariants

| Aspect | Behaviour |
|--------|-----------|
| `minLength` | minimum ORF length in **nucleotides**, inclusive (`Length ≥ minLength`); default 100 |
| `Position` | 0-based start offset within the scanned strand |
| `Frame` | reading-frame number 1–3 within the scanned strand |
| `IsReverseComplement` | true if found on the reverse complement (six-frame search) |
| `Sequence` / `Length` / `CodonCount` | span **including** the stop codon |
| null sequence | `ArgumentNullException`; shorter-than-a-codon / no-ATG / no-in-frame-stop → no ORFs |
| Case / alphabet | normalized by `DnaSequence`; non-ATG/non-stop codons are simply not boundary sites |

Invariants (INV-01..05): every `Sequence` starts `ATG`, ends TAA/TAG/TGA, `Length`
divisible by 3, `Length ≥ minLength`, `Frame ∈ {1,2,3}` with `IsReverseComplement`
selecting the strand. Complexity O(n²) worst case (many ATGs before a stop), O(n) typical.

## Worked oracles

- **Rosalind sample** `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`
  → the **4 distinct** proteins `MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`,
  `MTPRLGLESLLE` (last two share a stop).
- **Single forward ORF** `ATGAAAAAATAA` (minLength 1) → one ORF, `Sequence` `ATGAAAAAATAA`,
  position 0, frame 1, protein candidate `MKK`.

## Scope and the sibling annotation ORF finder

ATG-only, standard-code only; models neither codon bias, ribosome-binding sites, nor
splicing — ORF presence is not evidence of a real gene. It is deliberately **not
contract-equivalent** to the annotation-layer `GenomeAnnotator.FindOrfs` (test unit
ANNOT-ORF-001, primary spec `docs/algorithms/Annotation/ORF_Detection.md`), which
recognizes the prokaryotic start set ATG/GTG/TTG, measures `minLength` in **amino acids**
(default 100, excluding the terminal stop), exposes `searchBothStrands` (default true) /
`requireStartCodon` (default true) flags, emits 0-based half-open `[Start, End)` spans
whose `ProteinSequence` retains the terminal `*`, and accumulates **pending starts per
frame** so multiple starts sharing one stop each yield an ORF. Its companion
`GenomeAnnotator.FindLongestOrfsPerFrame` (internally calling `FindOrfs` with
`minLength: 1`) groups results under **signed frame keys** — `1..3` forward, `-1..-3`
reverse-complement — and keeps the longest translated ORF per frame. That annotation API
hard-codes `GeneticCode.Standard`, so genetic-code-parameterized ORF finding uses
`Translator.FindOrfs` instead. That annotation `FindOrfs` unit was independently
validated as [[annot-orf-001-report]] (Stage A PASS / Stage B PASS-WITH-NOTES, End state CLEAN,
35/35 tests, no defect — sole note is the non-canonical `requireStartCodon=false` run-off path). That annotation layer's ORF-based gene
prediction + Shine-Dalgarno RBS finder is [[prokaryotic-gene-prediction-rbs]]
(test unit ANNOT-GENE-001). Callers pick the entry point deliberately.
Coding-potential scoring of a candidate ORF is a separate step ([[coding-potential-hexamer-score]]).

## Reference sources

**Rosalind "ORF"** (definition, six-frame, distinct-protein return, nested-stop worked
example), **Wikipedia "Open reading frame"** (start/stop spans, six frame translations,
length-divisible-by-three), **NCBI ORFfinder** ("ATG only" default, nucleotide minimal-length
options, inclusive filter), and **NCBI Genetic Codes transl_table=1** (standard ATG start /
TAA-TAG-TGA stops). **No deviations** (one pre-existing greedy bug fixed); three source-anchored
assumptions — stop-inclusive span, nucleotide `minLength`, ATG-only start.
