---
type: concept
title: "Genetic code (codon → amino-acid translation tables)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/TRANS-CODON-001-Evidence.md
source_commit: 32cc779c6457154b8a1522281dacda677b854593
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
assigns. Validated as test unit **TRANS-CODON-001** ([[trans-codon-001-evidence]]); see
[[test-unit-registry]] for how the unit is tracked and [[algorithm-validation-evidence]] for the
evidence-artifact pattern.

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

## Scope

Four tables only (1/2/3/11) of NCBI's 33; single-codon lookup — whole-sequence translation with
frame handling lives in `Translator` / the six-frame ORF scanners
([[open-reading-frame-detection]]), and the MCP surface exposes it as `TranslateDna`/`TranslateRna`.

## Reference sources

**NCBI Genetic Codes** (`transl_table` 1/2/3/11 `AAs`+`Starts` strings, official spec),
**Wikipedia "Genetic code"** (64→20+3 degeneracy), **Wikipedia "Start codon"** (AUG universal,
GUG/UUG alternatives), **Wikipedia "Stop codon"** (UAA ochre / UAG amber / UGA opal), plus the
historical Nirenberg & Matthaei 1961 and Crick 1968 citations. Full record in
[[trans-codon-001-evidence]].
