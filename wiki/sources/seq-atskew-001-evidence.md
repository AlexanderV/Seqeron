---
type: source
title: "Evidence: SEQ-ATSKEW-001 (AT skew — (A−T)/(A+T) strand composition statistic)"
tags: [validation, sequence-statistics, composition]
doc_path: docs/Evidence/SEQ-ATSKEW-001-Evidence.md
sources:
  - docs/Evidence/SEQ-ATSKEW-001-Evidence.md
source_commit: ee3cfca8f1c41de229969aa234c2558284581909
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-ATSKEW-001

The validation-evidence artifact for test unit **SEQ-ATSKEW-001** — **AT skew**,
`(A − T) / (A + T)`, the A/T-strand-asymmetry sibling of GC skew. It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
[[test-unit-registry]] tracks the unit.

**No dedicated new concept for AT skew alone — it is one member of a family.** The formula,
range, zero-denominator rule, and biological meaning are synthesized on the family concept
[[nucleotide-composition-skew]] (which covers both AT skew and its GC-skew sibling); this
page records only what the artifact adds.

## What this file records

- **Online sources:**
  - **Lobry (1996)** *Mol Biol Evol* 13(5):660–665 (PMID 8676740; rank 1, primary) — the
    founding observation of intra-strand base asymmetry (*"departure from intrastrand
    equifrequency between A and T or between C and G"*) that skew quantifies; analyzed
    *E. coli*, *B. subtilis*, *H. influenzae*.
  - **Charneski et al. (2011)** *PLoS Genet* 7(9):e1002283 (rank 1) — the verbatim formula
    **AT skew = (A − T)/(A + T)**; Firmicute AT skew results from **selection, not
    mutation** (so AT skew need not track GC skew).
  - **Wikipedia "GC skew"** (rank 4, citing Lobry 1996) — corroborates both formulas and the
    value range: the composition-skew spectrum ranges **−1 … +1** (AT skew `−1 ⇔ A=0`,
    `+1 ⇔ T=0`).
  - **Biopython `Bio.SeqUtils.GC_skew`** (rank 3, reference implementation) — the
    symbol-handling conventions applied to the AT-skew analog: **case-insensitive** base
    counting (`count("G")+count("g")`), **zero-denominator ⇒ 0.0** (`ZeroDivisionError`
    caught), and **ambiguous/non-canonical bases ignored** (*"Does NOT look at any ambiguous
    nucleotides"*).
- **Datasets (hand-derived from the formula, arithmetic — no library run needed):**
  `AAAA → 1.0`, `TTTT → −1.0`, `ATAT → 0.0`, `AAAT → 0.5`, `ATTT → −0.5`,
  `GGCC → 0.0` (no A/T, zero denominator), `AAATGGGCCC → 0.5` (G/C ignored),
  `aaat → 0.5` (case-insensitive).
- **Corner cases / failure modes:** empty window / no A or T ⇒ `A+T = 0` ⇒ **0.0** (not
  NaN/exception); non-A/T symbols (gaps, `N`, IUPAC ambiguity) contribute to neither
  numerator nor denominator; bounds `[−1, +1]` with the two saturating endpoints above.

## Deviations and assumptions

**One documented assumption — lowercase + non-ACGT handling for the AT-skew member.**
Biopython ships a `GC_skew` line but **not** an AT-skew-specific one, so the
case-insensitive counting and "ignore everything that is not A/T" behaviour are taken by
analogy from the directly analogous `GC_skew` reference implementation rather than from an
AT-skew-specific source. The **formula itself is fully sourced** (Charneski 2011, Lobry
1996); only the symbol-handling convention is inferred by analogy, and it matches the
repository implementation (uppercases input via `ToUpperInvariant`, counts only `A`/`T`).

Recommended coverage (from the artifact): MUST — pure-A `+1.0` / pure-T `−1.0` (bounds);
balanced `A=T ⇒ 0.0` and asymmetric exact fractions (`AAAT ⇒ 0.5`); no A and no T `⇒ 0.0`
(zero denominator); G/C and other non-A/T symbols ignored. SHOULD — lowercase equals
uppercase; null string `⇒ 0` / empty `⇒ 0` / null `DnaSequence ⇒ ArgumentNullException`.
COULD — `DnaSequence` overload agrees with the string overload. No source contradictions —
Lobry, Charneski, Wikipedia, and the Biopython convention agree on formula and range.
