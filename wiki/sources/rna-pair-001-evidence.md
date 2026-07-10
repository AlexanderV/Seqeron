---
type: source
title: "Evidence: RNA-PAIR-001 (RNA base pairing — CanPair / GetBasePairType / GetComplement)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-PAIR-001-Evidence.md
sources:
  - docs/Evidence/RNA-PAIR-001-Evidence.md
source_commit: 59665e71035d6a10504325a61b2d3329858ebf36
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-PAIR-001

The validation-evidence artifact for test unit **RNA-PAIR-001** — **RNA Base Pairing**
(`RnaSecondaryStructure.CanPair` / `GetBasePairType` / `GetComplement`). It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
[[test-unit-registry]] tracks the unit.

This unit is the **RNA-secondary-structure family's own base-pairing primitive**: the pairing
predicate plus a **typed pair classifier** `GetBasePairType` (returns `WatsonCrick` / `Wobble` /
`null`) and the RNA base complement. It is the **same {A-U, G-C} + G-U wobble rule** already
synthesized on [[rna-base-pairing]] (the MiRNA family's `MiRnaAnalyzer` sibling, MIRNA-PAIR-001) —
so it **enriches that shared concept** rather than defining a new one. The classifier is the pairing
chemistry that the fold-scoring [[rna-free-energy-turner-model]] and the notation layer
[[rna-dot-bracket-notation]] both build on.

## What this file records

- **Canonical Watson-Crick pairs (RNA):** A•U (2 H-bonds) and G•C (3 H-bonds); in RNA thymine is
  replaced by uracil. `GetBasePairType` reports these as `WatsonCrick`.
- **G-U is a distinct Wobble pair, not Watson-Crick** (Crick 1966 wobble hypothesis): over the
  standard RNA alphabet, G pairs with C (WC) and U (wobble); U pairs with A (WC) and G (wobble).
  `GetBasePairType('G','U')` must report `Wobble`, never `WatsonCrick`. G-U is the *only* standard
  wobble over {A,C,G,U} (the inosine I-pairs are outside the standard alphabet).
- **Non-pairs → false / null:** any combination other than A-U, U-A, G-C, C-G, G-U, U-G does not pair
  (A-A, A-G, A-C, C-U, G-G, C-C → `CanPair` false, `GetBasePairType` null).
- **RNA complement (`GetComplement`):** A→U, U→A, G→C, C→G, and **T→A** (DNA T treated as U, whose
  complement is A); degenerate IUPAC complements preserved (N→N, R→Y, Y→R). Pairing itself is defined
  over the RNA alphabet {A,C,G,U}; T is a complement input, not a pairing input the sources define.
- **Sources:** Crick FHC (1966) *J Mol Biol* 19:548-555 (DOI 10.1016/S0022-2836(66)80022-0, the
  wobble hypothesis) + Wikipedia *Base pair* / *Wobble base pair* + IUPAC-IUB Commission (1970)
  *Biochemistry* 9(20):4022 (complement table) + Biopython `Bio.Seq.complement_rna` (reference impl,
  `complement_rna("CGAUT")` → `"GCUAA"`).

## Datasets / oracles

Canonical RNA base-pair truth table (Crick 1966; Wikipedia Base pair / Wobble base pair):

| base1 | base2 | CanPair | GetBasePairType |
|-------|-------|---------|-----------------|
| A | U | true | WatsonCrick |
| G | C | true | WatsonCrick |
| G | U | true | Wobble |
| A | A | false | null |
| C | U | false | null |

Complement (IUPAC / Biopython): A→U, U→A, G→C, C→G, T→A, N→N, R→Y, Y→R.

## Corner cases and invariants

- **Symmetry / order independence:** `CanPair(x,y) == CanPair(y,x)` and
  `GetBasePairType(x,y) == GetBasePairType(y,x)` for all pairs (pairing is reciprocal).
- **Wobble ⊄ Watson-Crick:** G-U is a real, thermodynamically comparable pair but a *separate* type.
- **Non-alphabet / out-of-range chars** return false / null with no exception (robustness, COULD-test).

## Deviations and assumptions

One recorded item — a **non-correctness-affecting normalization**: the implementation upper-cases
inputs before lookup (case-insensitive pairing); the sources define pairing over uppercase bases, and
case denotes the same nucleotide, so this does not change which base is meant. All pairing-type and
complement values are source-backed. **No source contradictions** — Crick 1966, Wikipedia,
IUPAC-IUB 1970, and Biopython agree.
