---
type: source
title: "Evidence: TRANS-CODON-001 (Codon translation, genetic-code tables 1/2/3/11)"
tags: [validation, annotation]
doc_path: docs/Evidence/TRANS-CODON-001-Evidence.md
sources:
  - docs/Evidence/TRANS-CODON-001-Evidence.md
source_commit: 32cc779c6457154b8a1522281dacda677b854593
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-CODON-001

The validation-evidence artifact for test unit **TRANS-CODON-001** — **codon → amino-acid
translation** via the NCBI genetic-code tables (`GeneticCode.Translate`, area *Translation*).
It is one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; the algorithm, its four supported tables, contract, degeneracy, and corner
cases are synthesized in [[genetic-code-translation]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources:**
  - **NCBI Genetic Codes** (`wprintgc.cgi`, official reference) — complete `transl_table`
    1–33; the artifact reproduces the `AAs` and `Starts` strings for **tables 1, 2, 3, 11**
    (updated 2024-09-23) as the exact oracle.
  - **Wikipedia "Genetic code"** — 64 codons → 20 amino acids + 3 stops, degeneracy.
  - **Wikipedia "Start codon"** — AUG universal start; GUG/UUG alternatives in bacteria.
  - **Wikipedia "Stop codon"** — UAA (ochre), UAG (amber), UGA (opal) nomenclature.
  - Historical: Nirenberg & Matthaei 1961; Crick 1968.
- **Datasets (documented oracles):**
  - **Complete Standard table** — all 64 codons → single-letter AA / `*` (the primary
    coverage fixture).
  - **Table-2 (Vertebrate Mito) diffs**: AGA→`*`, AGG→`*`, AUA→`M`, UGA→`W` (vs Standard
    R/R/I/`*`).
  - **Table-3 (Yeast Mito) diffs**: CUU/CUC/CUA/CUG→`T`, AUA→`M`, UGA→`W`.
  - Degeneracy table (Met/Trp single-codon; Leu/Ser/Arg six-fold; 3 stops).
- **Corner cases / failure modes documented**: codon length ≠ 3 → exception; DNA↔RNA both
  work (T↔U); case-insensitive; stop → `'*'`; start-codon-mid-sequence → `M` (not fMet);
  empty/short/long/null/invalid-nucleotide/invalid-table → `ArgumentException` /
  `ArgumentNullException`.

## Deviations and assumptions

The Evidence doc's own *Deviations and Assumptions* section reads **"None. Implementation
matches NCBI translation tables exactly"** — accurate for the **codon→AA mappings and
start/stop sets** (derived verbatim from the NCBI `AAs`/`Starts` strings for tables 1/2/3/11).

**One source-vs-implementation discrepancy flagged** (API-contract layer, not the code tables):
the doc's *Documented Corner Cases* and *Known Failure Modes* say an "Unknown codon (e.g.,
NNN)" / "Invalid nucleotide (X, Z)" should raise `ArgumentException`. The actual
`GeneticCode.Translate` (`GeneticCode.cs`) returns **`'X'`** for any **valid IUPAC ambiguity
codon** (alphabet `ACGURYMKSWBDHVN`, so `NNN`, `RAY`, … are untranslatable-but-not-invalid) and
throws only for a triplet containing a **non-IUPAC** character (e.g. `Z`). So `NNN → 'X'`, not
an exception. Recorded on [[genetic-code-translation]] for reconciliation.

## Implementation note

Seqeron supports **4** of NCBI's 33 tables: Standard (1), Vertebrate Mitochondrial (2), Yeast
Mitochondrial (3), Bacterial/Archaeal/Plant-Plastid (11), as static singletons
(`GeneticCode.Standard` etc.) plus `GetByTableNumber`. Alongside the single-codon `Translate`,
the surface exposes `IsStartCodon`, `IsStopCodon`, and the reverse `GetCodonsForAminoAcid`.
This is the genetic-code table that [[open-reading-frame-detection]] and [[codon-optimization]]
read; whole-sequence framed translation lives in `Translator`.
