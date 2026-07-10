---
type: source
title: "Evidence: SEQ-RNACOMP-001 (RNA-specific per-base complement — GetRnaComplementBase, IUPAC-complete)"
tags: [validation, rna]
doc_path: docs/Evidence/SEQ-RNACOMP-001-Evidence.md
sources:
  - docs/Evidence/SEQ-RNACOMP-001-Evidence.md
source_commit: 51ed4d23872ce7c6646683d002e13e9388412d53
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-RNACOMP-001

The validation-evidence artifact for test unit **SEQ-RNACOMP-001** — **RNA-specific Complement**
(`GetRnaComplementBase`), the **per-base, IUPAC-complete** RNA complement lookup. It is one instance
of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
[[test-unit-registry]] tracks the unit.

> **Naming caveat (RNACOMP ≠ composition):** despite the `RNACOMP` slug this unit is **RNA
> complement**, *not* RNA base composition. It has nothing to do with the composition tally in
> [[base-composition]]; its subject is the single-base complement map, so it enriches
> [[rna-base-pairing]] (which already owns the RNA base complement).

This is the RNA sibling of the DNA per-base complement `GetComplementBase` (**SEQ-COMP-001**, not yet
ingested). It is the **full-IUPAC-complete** per-base complement over the RNA alphabet — the same base
complement chemistry as `RnaSecondaryStructure.GetComplement` ([[rna-base-pairing]], RNA-PAIR-001) but
a distinct SEQ-\* sequence-utility surface that maps **all** IUPAC ambiguity codes, not just the
canonical bases.

## What this file records

- **Canonical RNA complement:** A→U, C→G, G→C, U→A. In RNA, thymine is replaced by uracil.
- **T treated as U:** DNA `T` in an RNA context complements to **A** (Biopython builds the RNA table
  as `ambiguous_rna_complement["T"] = ...["U"] = "A"`; docstring: "Any T in the sequence is treated as
  a U"). So both `T` and `U` → A.
- **Full IUPAC ambiguity map (RNA alphabet):** R↔Y, M↔K, D↔H, B↔V reciprocal pairs; **self-complementary**
  W→W, S→S, X→X, N→N. Every code's complement is the code of the complemented base set
  (e.g. R=A|G → Y=C|U). The RNA table is **identical to the DNA `ambiguous_dna_complement` table except
  the base alphabet swaps T→U**.
- **Repo casing convention:** `GetRnaComplementBase` returns **uppercase** for recognized bases
  (the DnaSequence/RnaSequence normalize-to-uppercase convention, shared with SEQ-COMP-001). This is the
  **only** behavioral divergence from Biopython, which preserves input case (`a → u`). The divergence is
  casing-only — it does not change which base pairs with which.
- **Pass-through of non-IUPAC characters:** letters/symbols not in the table (e.g. `Z`, gaps `-`/`.`,
  digits) are **returned unchanged**, not an error; gap → gap.
- **Involution:** complement is an involution on the canonical RNA bases and ambiguity codes (within the
  U-alphabet; `T` is absorbed into U so `complement(complement(T)) = complement(A) = U`, not T).

## Datasets / oracles

Full-alphabet worked example (Biopython 1.79 `reverse_complement_rna("ACGTUacgtuXYZxyz")` →
`'zrxZRXaacguAACGU'`; un-reversed forward complement `"UGCAAugcaaXRZxrz"`), per input char with the
repo's uppercase convention: A→U, C→G, G→C, T→A, U→A, X→X (pass-through), Y→R, Z→Z (pass-through);
lowercase recognized bases uppercase to the same complement (a→U, g→C, y→R), lowercase non-bases pass
through as-is (x→x, z→z).

IUPAC ambiguity-code complements (RNA): R(A,G)→Y, Y(C,U)→R, S(G,C)→S, W(A,U)→W, K(G,U)→M, M(A,C)→K,
B(C,G,U)→V, D(A,G,U)→H, H(A,C,U)→D, V(A,C,G)→B, N(any)→N.

## Corner cases and invariants

- **RNA-vs-DNA distinction (COULD-test):** `GetRnaComplementBase('A') = 'U'` vs
  `GetComplementBase('A') = 'T'` — confirms the RNA path is not the DNA path.
- **T → A** (not preserved as T) is the documented RNA-context behavior.
- **Symmetry of reciprocal codes:** A↔U, C↔G, R↔Y, M↔K, D↔H, B↔V; self-complementary W/S/X/N.

## Deviations and assumptions

One recorded item — a **non-correctness-affecting normalization**: recognized bases are returned
**uppercase** (repo convention, per SEQ-COMP-001 MUST-02 and the method's XML remarks), whereas
Biopython preserves case. Casing only; the complement identity is unchanged. All complement values are
source-backed. **No source contradictions** — Biopython (`IUPACData.py`, `Seq.py`, docs),
bioinformatics.org SMS, and the NC-IUB 1984 standard (Cornish-Bowden 1985, *NAR* 13(9):3021) agree.
