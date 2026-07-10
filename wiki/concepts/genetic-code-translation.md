---
type: concept
title: "Genetic code (codon → amino-acid translation tables)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/TRANS-CODON-001-Evidence.md
  - docs/Evidence/TRANS-PROT-001-Evidence.md
source_commit: 7122c87a12f8c52d64ed7d5f5241ff5aa19879ef
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-codon-001-evidence
      evidence: "Test Unit ID: TRANS-CODON-001 — Area: Translation — codon → amino-acid lookup (GeneticCode.Translate)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-prot-001-evidence
      evidence: "Test Unit ID: TRANS-PROT-001 — Area: Translation — whole-sequence framed / six-frame translation + ORF (Translator)"
      confidence: high
      status: current
---

# Genetic code (codon → amino-acid translation tables)

The **genetic code** is the many-to-one map from **64 nucleotide triplets (codons)** to
**20 amino acids + 3 stop signals**. Seqeron models it as `GeneticCode` in
`Seqeron.Genomics.Core` (`GeneticCode.cs`): a single-codon lookup (`Translate(codon) → char`)
plus start/stop predicates and a reverse lookup. It is the **foundational table** that every
higher codon-level operation reads — [[open-reading-frame-detection|ORF detection]] scans it
for ATG→stop spans, [[codon-optimization]] substitutes within its synonymous families, and the
[[relative-synonymous-codon-usage|RSCU]] / [[codon-adaptation-index|CAI]] /
[[effective-number-of-codons|ENC]] measures all partition codons by the amino acid this table
assigns. Validated as test unit **TRANS-CODON-001** ([[trans-codon-001-evidence]]); the
**whole-sequence** layer above it — framed / six-frame translation and ORF finding via the
`Translator` class — is validated as **TRANS-PROT-001** ([[trans-prot-001-evidence]]) and
described below. See [[test-unit-registry]] for how the units are tracked and
[[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The four supported tables

Seqeron ships **four** of the NCBI genetic-code tables (out of 33), each a static singleton
(`GeneticCode.Standard`, `.VertebrateMitochondrial`, `.YeastMitochondrial`, `.BacterialPlastid`,
plus `GetByTableNumber(int)`):

| NCBI table | Property | Key deviations from Standard |
|-----------|----------|------------------------------|
| **1** Standard | `Standard` | baseline; starts AUG (primary) + UUG/CUG (alt) |
| **2** Vertebrate Mitochondrial | `VertebrateMitochondrial` | AGA/AGG → **Stop** (not Arg); AUA → **Met** (not Ile); UGA → **Trp** (not Stop) |
| **3** Yeast Mitochondrial | `YeastMitochondrial` | CUN (CUU/CUC/CUA/CUG) → **Thr** (not Leu); AUA → **Met**; UGA → **Trp** |
| **11** Bacterial/Archaeal/Plant-Plastid | `BacterialPlastid` | same codon table as Standard, **extra start codons** (GUG/UUG/CUG/AUU/AUC/AUA) |

Codon→amino-acid mappings and the start/stop sets are taken **directly from the NCBI `AAs` and
`Starts` strings** for tables 1/2/3/11 (NCBI Genetic Codes, updated 2024-09-23); the Evidence
records **no deviation** in the mappings themselves. Stop codons translate to `'*'`.

## Degeneracy (codon redundancy)

The code is degenerate: most amino acids have several synonymous codons. Only **Met (AUG)** and
**Trp (UGG)** are single-codon families — the same fact that makes them fixed points under
[[codon-optimization]] and gives them relative adaptiveness `w ≡ 1` in [[codon-adaptation-index|CAI]].
Six-fold families (Leu, Ser, Arg) sit at the other extreme. `GetCodonsForAminoAcid(aa)` is the
reverse lookup used to verify degeneracy.

## Input contract and normalization

- **DNA or RNA** input both accepted — `Translate` normalizes `T → U` internally (so `ATG`≡`AUG`).
- **Case-insensitive** (`AUG` = `aug` = `AuG`, upper-cased internally).
- Codon **must be exactly 3 characters**; empty/`null`/wrong-length → `ArgumentException`
  (`ArgumentNullException` for null).
- `IsStartCodon` / `IsStopCodon` test membership in the table's start/stop sets; a start codon
  used mid-sequence still translates to its amino acid (**M**, not fMet — the code does not model
  formyl-methionine).

## IUPAC ambiguity → 'X' (a deviation from the Evidence doc's corner-case table)

The implementation treats **valid IUPAC ambiguity codons** (any triplet over
`ACGURYMKSWBDHVN`, e.g. `NNN`, `RAY`) as **untranslatable but not invalid**: `Translate`
returns `'X'` (unknown amino acid) rather than throwing. Only a triplet containing a
**non-IUPAC** character (e.g. `Z`) raises `ArgumentException`. This **contradicts** the
Evidence doc's *Documented Corner Cases* / *Known Failure Modes* tables, which state that an
"Unknown codon (e.g., NNN)" and "Invalid nucleotide (X, Z)" should both yield an
`ArgumentException`. The mapping tables match NCBI exactly; this divergence is in the
**API-contract layer** (ambiguity handling), not in the code tables. Flagged as a
source-vs-implementation discrepancy for reconciliation. (`GeneticCode.cs`, `Translate`.)

## Whole-sequence framed translation (the `Translator` layer)

Above the single-codon lookup sits the **`Translator`** class (test unit **TRANS-PROT-001**,
[[trans-prot-001-evidence]]), which walks a whole nucleotide sequence in triplets and emits a
protein. It **composes** this genetic-code table rather than duplicating it — the code table is
selectable via a `GeneticCode` parameter, so all four tables (1/2/3/11) and their alternative
start/stop sets flow through unchanged.

- **Reading frame.** A `frame` parameter **0, 1, or 2** offsets the first codon; other values
  throw. Each single strand thus has three frames, and reverse-complementing the strand gives
  three more — **six frames** total (Wikipedia "Reading frame").
- **Six-frame translation** returns a dictionary keyed **−3…−1, +1…+3 (0 excluded)**: three
  forward frames and three on the reverse complement, each read 5'→3'.
- **`toFirstStop`** optionally terminates translation at the first in-frame stop (the codon-table
  `'*'`), matching the "translate until a stop" convention; otherwise stops render as `*` and
  reading continues to the sequence end.
- **Input contract** mirrors the codon layer: **DNA or RNA** (automatic `T→U`), case-insensitive;
  a trailing partial codon (length not divisible by 3) is simply not translated.
- **ORF finding.** `Translator.FindOrfs` is a **genetic-code-parameterized** ORF scanner
  (configurable minimum length, both-strand search) — deliberately **not** contract-equivalent to
  the ATG-only / standard-code `GenomicAnalyzer.FindOpenReadingFrames`
  ([[open-reading-frame-detection]]) nor the prokaryotic-start annotation-layer
  `GenomeAnnotator.FindOrfs`; callers pick the entry point deliberately.

**Correctness oracle:** the human insulin **B chain** (UniProt P01308, positions 25–54) DNA
`TTCGTG…AAGACC` (90 nt) → `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` (30 aa). All four tables verified
codon-by-codon against NCBI (2024-09-23); **no deviation** recorded for translation itself.

## Scope

Four tables only (1/2/3/11) of NCBI's 33. The single-codon lookup is `GeneticCode.Translate`;
whole-sequence framed / six-frame translation and ORF finding live in `Translator` (above), and
the MCP surface exposes translation as `TranslateDna`/`TranslateRna`.

## Reference sources

**NCBI Genetic Codes** (`transl_table` 1/2/3/11 `AAs`+`Starts` strings, official spec),
**Wikipedia "Genetic code"** (64→20+3 degeneracy), **Wikipedia "Start codon"** (AUG universal,
GUG/UUG alternatives), **Wikipedia "Stop codon"** (UAA ochre / UAG amber / UGA opal), plus the
historical Nirenberg & Matthaei 1961 and Crick 1968 citations. Full record in
[[trans-codon-001-evidence]]. For the whole-sequence layer: **Wikipedia "Translation (biology)"**
(triplet reading, start/stop, release factors), **Wikipedia "Reading frame"** (three/six frames),
**Wikipedia "Open reading frame"** (six-frame, sORFs), and **UniProt P01308** (insulin B-chain
oracle) — full record in [[trans-prot-001-evidence]].
